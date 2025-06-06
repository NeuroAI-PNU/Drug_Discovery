import time
import torch
import shutil
import pathlib
import subprocess
import numpy as np
from argparse import Namespace
import pandas as pd
from datetime import datetime
from os.path import join
from tqdm import tqdm
from pathlib import Path
from torch.optim import Adam, SGD, lr_scheduler
from ranger import Ranger
from torch_geometric.data import DataLoader
from model import Model
from metrics import regression_metrics, binary_metrics_multi_target_nan
from metrics import blend_regression, blend_binary_classification_mt, auto_metrics
from loss import get_loss
from model import model_args
from dataset import auto_dataset
from utils import GPUManager
from logger import auto_summarize_logs, config2cmd, read_logs
from dataset import perturb_test, dataset_names

torch.backends.cudnn.enabled = True


class Trainer():
    def __init__(self, args, model, train_dataset, valid_dataset, test_dataset=None, print_log=True):
        self.args = args  # ;  print(args)
        self.device = torch.device("cuda:{}".format(
            args.gpu) if torch.cuda.is_available() else "cpu")
        self.model = model.to(self.device)

        # Setting the train valid and test data loader
        self.train_dataloader = DataLoader(
            train_dataset, batch_size=args.batch_size)
        self.valid_dataloader = DataLoader(valid_dataset, batch_size=32)
        if test_dataset:
            self.test_dataloader = DataLoader(test_dataset, batch_size=32)

        # Setting the loss, optimizer, scheduler
        self.criterion = get_loss(args.loss)
        if args.optim == 'Ranger':
            # weight_decay=args.weight_decay
            self.optimizer = Ranger(
                self.model.parameters(), lr=args.lr, k=args.k)
        elif args.optim == 'Adam':
            self.optimizer = Adam(self.model.parameters(), lr=args.lr)
        elif args.optim == 'SGD':
            self.optimizer = SGD(self.model.parameters(), lr=args.lr)
        else:
            raise Exception('Error optimizer argv')
        self.scheduler = lr_scheduler.ReduceLROnPlateau(self.optimizer, mode='min', factor=args.lr_reduce_rate,
                                                        patience=args.lr_reduce_patience, min_lr=1e-6)

        # other
        self.start = time.time()
        self.print_log = print_log
        cur_abs_file_dir = Path(__file__).resolve().parent
        save_id = datetime.utcnow().strftime(
            '%Y-%m-%d_%H:%M:%S.%f')[:-3] + '_' + 'seed_' + str(args.seed)
        self.log_save_dir = cur_abs_file_dir / \
            ('log_' + self.args.dataset) / '{}'.format(save_id)
        Path.mkdir(self.log_save_dir.parent, exist_ok=True)
        Path.mkdir(self.log_save_dir, exist_ok=True)
        self.records = {'val_losses': []}

        self.log(msgs=['\t{}:{}\n'.format(k, v)
                 for k, v in args.__dict__.items()])
        self.log('save id: {}'.format(save_id))
        self.log('run device: {}'.format(self.device))
        self.log('train set num:{}    valid set num:{}    test set num: {}'.format(
            len(train_dataset), len(valid_dataset), len(test_dataset)))
        self.log("total parameters:" + str(sum([p.nelement() for p in self.model.parameters()])))
        self.log(msgs=str(model).split('\n'))

    def train(self):
        self.log('Training start...')
        early_stop_cnt = 0
        for epoch in tqdm(range(self.args.epochs)):
            trn_loss = self.train_iterations()
            val_loss, result = self.valid_iterations()
            self.scheduler.step(val_loss)
            lr_cur = self.scheduler.optimizer.param_groups[0]['lr']
            self.log('Epoch:{} trn_loss:{:.5f} val_loss:{:.5f} val_result:{} lr_cur:{:.7f}'.format(
                epoch, trn_loss, val_loss, result, lr_cur), with_time=True)
            self.records['val_losses'].append(val_loss)
            if val_loss == np.array(self.records['val_losses']).min():
                self.save_ckpt(epoch)
                early_stop_cnt = 0
            else:
                early_stop_cnt += 1
            if 0 < self.args.early_stop_patience < early_stop_cnt:
                self.log('Early stop hitted!')
                break
        self.save_ckpt(epoch, final_save=True)

    def train_and_test(self):
        self.train()
        self.log('Testing...')
        self.load_best_ckpt()
        val_loss, val_result = self.valid_iterations(mode='valid')
        test_loss, test_result = self.valid_iterations(mode='test')
        self.log(msg=str(vars(self.args)))
        loss_info = {'testloss': test_loss.item(), 'valloss': val_loss.item()}
        val_result_new = {}
        for k in val_result.keys():
            val_result_new['val' + k] = val_result[k]
        self.log('{}|{}|{}'.format(loss_info, test_result, val_result_new))

    def save_ckpt(self, epoch, final_save=False):
        file_name = 'final_save.ckpt' if final_save else 'best_save.ckpt'
        with open(join(self.log_save_dir, file_name), 'wb') as f:
            torch.save({
                'args': self.args,
                'records': self.records,
                'model_state_dict': self.model.state_dict(),
            }, f)
        self.log('Model saved at epoch {}'.format(epoch))

    def gen_test_batch(self):
        for batch in self.valid_dataloader:
            torch.save(batch, './other/test_batch.pt')
            break

    def load_best_ckpt(self):
        ckpt_path = join(self.log_save_dir, 'best_save.ckpt')
        self.log('The best ckpt is {}'.format(ckpt_path))
        self.load_ckpt(ckpt_path)

    def load_ckpt(self, ckpt_path):
        self.log('Ckpt loading: {}'.format(ckpt_path))
        ckpt = torch.load(ckpt_path, map_location=self.device)
        self.args = ckpt['args']
        self.records = ckpt['records']
        self.model.load_state_dict(ckpt['model_state_dict'])

    def log(self, msg=None, msgs=None, with_time=False):
        if self.print_log is False:
            return
        if with_time:
            msg = msg + ' time elapsed {:.2f} hrs ({:.1f} mins)'.format(
                (time.time() - self.start) / 3600.,
                (time.time() - self.start) / 60.
            )
        with open(join(self.log_save_dir, 'log.txt'), 'a+') as f:
            if msgs:
                self.log('#' * 80)
                if '\n' not in msgs[0]:
                    msgs = [m + '\n' for m in msgs]
                f.writelines(msgs)
                for x in msgs:
                    print(x, end='')
                self.log('#' * 80)
            if msg:
                f.write(msg + '\n')
                print(msg)


class TrainerMolBinaryClassificationNAN(Trainer):
    def __init__(self, args, model, train_dataset, valid_dataset, test_dataset=None, print_log=True):
        super(TrainerMolBinaryClassificationNAN, self).__init__(args, model, train_dataset, valid_dataset, test_dataset,
                                                                print_log)
        self.metrics_fn = binary_metrics_multi_target_nan

    def train_iterations(self):
        self.model.train()
        losses = []
        for i, data in enumerate(self.train_dataloader):
            if data.__sizeof__() <= 1:
                continue  # bn will raise error if there is only one GLAM
            self.optimizer.zero_grad()
            mol_batch = data.to(self.device)
            y_true = data.y.to(self.device)
            y_scores = self.model(mol_batch)
            # (N, T*C) to (N, T, C)  where C is n_class
            y_scores = y_scores.view(y_true.shape[0], y_true.shape[1], 2)
            y_scores = y_scores.view([-1, 1])
            valid_id = torch.where(y_true[:, :] >= 0)
            #loss = self.criterion(y_scores[valid_id], y_true[valid_id])  # (1)
            loss = self.criterion(y_scores, y_true)  # (1)
            loss.backward()
            self.optimizer.step()
            losses.append(loss.item())
            if i % self.args.verbose_patience == 0:
                self.log('\tbatch {} training loss: {:.5f}'.format(
                    i, loss.item()), with_time=True)
        trn_loss = np.array(losses).mean()
        return trn_loss

    @torch.no_grad()
    def valid_iterations(self, mode='valid'):
        self.model.eval()
        dataloader = self.valid_dataloader if mode == 'valid' else self.test_dataloader
        losses, ys_true, ys_score, ys_pred = [], [], [], []
        for data in dataloader:
            if data.__sizeof__() <= 1:
                continue  # bn will raise error if there is only one GLAM
            mol_batch = data.to(self.device)
            y_true = data.y.to(self.device)
            y_scores = self.model(mol_batch)  # (N, T*C)
            # (N, T*C) to (N, T, C)  where C is n_class
            y_scores = y_scores.view(y_true.shape[0], y_true.shape[1], 2)
            valid_id = torch.where(y_true[:, :] >= 0)  # (N, T)
            # loss = self.criterion(y_scores[valid_id], y_true[valid_id])  # (1)
            loss = self.criterion(y_scores, y_true)  # (1)
            losses.append(loss.item())

            # (N, T, C) -> (N, T, C) -> (N, T)
            y_score = torch.softmax(y_scores, 2)[:, :, 1]
            # dim 2 is C class     # (N, T, C) -> (N, T)
            y_pred = torch.argmax(y_scores, dim=2)
            ys_pred.append(y_pred.cpu())
            ys_score.append(y_score.cpu())
            ys_true.append(y_true.cpu())
        mean_loss = torch.tensor(losses).mean()  # (N) -> (1)
        ys_true, ys_score, ys_pred = torch.cat(ys_true, dim=0), torch.cat(
            ys_score, dim=0), torch.cat(ys_pred, dim=0)
        if mode == 'inference':
            return ys_true, ys_score, ys_pred  # (N, T), (N, T), (N, T)
        result = self.metrics_fn(
            ys_true.numpy(), ys_score.numpy(), ys_pred.numpy())
        return mean_loss, result


class TrainerMolBinaryClassificationNANBCE(Trainer):
    def __init__(self, args, model, train_dataset, valid_dataset, test_dataset=None, print_log=True):
        super(TrainerMolBinaryClassificationNANBCE, self).__init__(args, model, train_dataset, valid_dataset,
                                                                   test_dataset,
                                                                   print_log)
        self.metrics_fn = binary_metrics_multi_target_nan

    def train_iterations(self):
        self.model.train()
        losses = []
        for i, data in enumerate(self.train_dataloader):
            if data.__sizeof__() <= 1:
                continue  # bn will raise error if there is only one GLAM
            self.optimizer.zero_grad()
            mol_batch = data.to(self.device)
            y_score = self.model(mol_batch)
            y_true = data.y.to(self.device)
            loss = self.criterion(
                y_score[y_true >= 0], y_true[y_true >= 0].float())
            loss.backward()
            self.optimizer.step()
            losses.append(loss.item())
            if i % self.args.verbose_patience == 0:
                self.log('\tbatch {} training loss: {:.5f}'.format(
                    i, loss.item()), with_time=True)
        trn_loss = np.array(losses).mean()
        return trn_loss

    @torch.no_grad()
    def valid_iterations(self, mode='valid'):
        self.model.eval()
        dataloader = self.valid_dataloader if mode == 'valid' else self.test_dataloader
        losses, ys_true, ys_score = [], [], []
        for data in dataloader:
            if data.__sizeof__() <= 1:
                continue  # bn will raise error if there is only one GLAM
            mol_batch = data.to(self.device)
            y_score = self.model(mol_batch)
            y_true = data.y.to(self.device)
            loss = self.criterion(
                y_score[y_true >= 0], y_true[y_true >= 0].float())
            losses.append(loss.item())
            ys_score.append(torch.sigmoid(y_score).cpu())  # add sigmoid
            ys_true.append(y_true.cpu())
        mean_loss = torch.tensor(losses).mean()
        ys_true, ys_score = torch.cat(
            ys_true, dim=0), torch.cat(ys_score, dim=0)
        if mode == 'inference':
            return ys_score, ys_true
        result = self.metrics_fn(ys_true.numpy(), ys_score.numpy())
        return mean_loss, result
    
    @torch.no_grad()
    def valid_iterations_only_predresult(self, mode='valid'):
        self.model.eval()
        dataloader = self.valid_dataloader if mode == 'valid' else self.test_dataloader
        ys_score = [] # losses = [] ,ys_true = [] were removed
        for data in dataloader:
            if data.__sizeof__() <= 1:
                continue  # bn will raise error if there is only one GLAM
            mol_batch = data.to(self.device)
            y_score = self.model(mol_batch)
            #y_true = data.y.to(self.device)
            #loss = self.criterion(
                #y_score[y_true >= 0], 
                #y_true[y_true >= 0].float())
            #losses.append(loss.item())
            ys_score.append(torch.sigmoid(y_score).cpu())  # add sigmoid
            #ys_true.append(y_true.cpu())
        #mean_loss = torch.tensor(losses).mean()
        #ys_true, ys_score = torch.cat(
            #ys_true, dim=0), torch.cat(ys_score, dim=0)
        if mode == 'inference':
            return ys_score #ys_true
        #result = self.metrics_fn(ys_score.numpy()) # ys_true.numpy() was removed
        #result = ys_score
        #return result #mean_loss was removed


class TrainerMolRegression(Trainer):
    def __init__(self, args, model, train_dataset, valid_dataset, test_dataset=None, print_log=True):
        super(TrainerMolRegression, self).__init__(args, model,
                                                   train_dataset, valid_dataset, test_dataset, print_log)
        self.metrics_fn = regression_metrics

    def train_iterations(self):
        self.model.train()
        losses = []
        for i, data in enumerate(self.train_dataloader):
            if data.__sizeof__() <= 1:
                continue  # bn will raise error if there is only one GLAM
            self.optimizer.zero_grad()
            mol_batch = data.to(self.device)
            y_true = data.y.to(self.device).view(-1)
            output = self.model(mol_batch).view(-1)
            loss = self.criterion(output, y_true)
            loss.backward()
            self.optimizer.step()
            losses.append(loss.item())
            if i % self.args.verbose_patience == 0:
                self.log('\tbatch {} training loss: {:.5f}'.format(
                    i, loss.item()), with_time=True)
        trn_loss = np.array(losses).mean()
        return trn_loss

    @torch.no_grad()
    def valid_iterations(self, mode='valid'):
        self.model.eval()
        dataloader = self.valid_dataloader if mode == 'valid' else self.test_dataloader
        losses, ys_true, ys_pred = [], [], []
        for data in dataloader:
            if data.__sizeof__() <= 1:
                continue  # bn will raise error if there is only one GLAM
            mol_batch = data.to(self.device)
            y_true = data.y.to(self.device).view(-1)
            output = self.model(mol_batch).view(-1)
            loss = self.criterion(output, y_true)
            losses.append(loss.cpu())
            ys_true.append(y_true.cpu())
            ys_pred.append(output.cpu())

        mean_loss = torch.tensor(losses).mean()
        if mode == 'inference':
            return torch.cat(ys_true), torch.cat(ys_pred)
        result = self.metrics_fn(
            torch.cat(ys_true).numpy(), torch.cat(ys_pred).numpy())
        return mean_loss, result

    def pasp(self):
        for level in [1, 2, 3]:
            self.log('Run model for perturbed test level {}...'.format(level))
            M, M_prime, Q, Q_prime = perturb_test(self.args.dataset, level)

            self.test_dataloader = DataLoader(M, batch_size=32)
            self.valid_iterations(mode='test')
            _, P = self.valid_iterations(mode='inference')
            self.test_dataloader = DataLoader(M_prime, batch_size=32)
            self.valid_iterations(mode='test')
            _2, P_prime = self.valid_iterations(mode='inference')

            L_P_P_prime = regression_metrics(P.numpy(), y_pred=P_prime.numpy())
            L_Q_Q_prime = regression_metrics(Q, y_pred=Q_prime)
            self.log('L(P, P\') is {}, and\n L(Q, Q\') is {}'.format(
                L_P_P_prime, L_Q_Q_prime))
            self.log('\Delta_RMSE={}'.format(
                L_P_P_prime['rmse'] - L_Q_Q_prime['rmse']))


def _blend_and_inference(ids, configs, custom_dataset, log=print):
    from dataset import Dataset
    log('Start to blend models and inference ...')
    #outputs = []
    for id, config in zip(ids, configs):
        args = Namespace(**eval(config))
        args.gpu = 0  # all
        args, dataset, _Trainer = auto_dataset(args)
        model = Model(dataset.mol_num_node_features,
                      dataset.mol_num_edge_features, **model_args(args))
        trainer = _Trainer(args, model, dataset.train,
                           dataset.val, dataset.test, print_log=False)
        # remove new made log directory
        shutil.rmtree(trainer.log_save_dir)
        # print(trainer.test_dataloader.dataset[0].smi)
        
        # 여기에 custom_data에 관한 코드를 짜서 넣어주자 tanimoto -> pk object로 만드는 코드
        if custom_dataset is not None:
            processed_custom_dataset = Dataset(args.dataset_root, dataset = custom_dataset, split_seed= args.split_seed)
            print(processed_custom_dataset)
            trainer.test_dataloader = DataLoader(
                processed_custom_dataset, batch_size=32)
            log('Customed test dataset loaded!')
        #print(trainer.test_dataloader.dataset[0].smi)
        trainer.log_save_dir = './log_{}/{}'.format(args.dataset, id)
        trainer.load_best_ckpt()
        log('Checkpoint {} loaded.'.format(id))
        output = trainer.valid_iterations_only_predresult(mode='inference') # was changed from valid_interations to valid_iterations_only_predresult
        #print(output)
        #print(type(output))
        smiles = []
        for i in range(len(trainer.test_dataloader.dataset)):
            smile = trainer.test_dataloader.dataset[i].smi
            smiles.append(smile)
        import pickle
        with open('./smiles.pk', 'wb') as f:
            pickle.dump(smiles, f)
        with open('./output_.pk', 'wb') as f:
            pickle.dump(output, f)
        #outputs.append(output)
        log('inference done!', with_time=True)
    #log('blend results: ')
    #if args.dataset in dataset_names["r"]:  # + ['physprop_perturb']:
    #    log(blend_regression(outputs))
    #elif args.dataset in dataset_names["c"] + ['demo']:
    #    log(blend_binary_classification_mt(
    #        outputs, metrics_fn=binary_metrics_multi_target_nan))
    #else:
        #raise ValueError('unknown dataset')
    #if args.dataset in ['physprop_perturb']:
        #return blend_regression(outputs, return_pred=True)
    log('Done!', with_time=True)


class GLAMHelper():
    def __init__(self, dataset: str, n_blend=3):
        self.dataset = dataset
        self.n_blend = n_blend
        self.start = time.time()
        self.logs_dir = pathlib.Path('./log_{}/'.format(self.dataset))
        self.log('GLAMHelper pr of {} start...'.format(dataset))

    def blend_and_inference(self, custom_dataset = None):
        ids, configs = self.select_top_config()
        custom_dataset = custom_dataset
        #print(custom_dataset)
        return _blend_and_inference(ids, configs, custom_dataset, log=self.log)

    def select_top_config(self):
        logs = read_logs(self.logs_dir)
        if len(logs) < 1:
            self.log(
                'Error: There is no log files found in {}!'.format(self.logs_dir))
        logs_pd = pd.DataFrame(logs)
        # index = logs_pd.loc[:, 'epochs'] >= 100  # assert evaluation epochs > 100 > search epochs
        metrics = auto_metrics(self.dataset)
        n_blend = len(logs_pd) if len(logs_pd) < self.n_blend else self.n_blend
        self.log('{} checkpoints select!'.format(n_blend))
        logs_pd_selected = logs_pd.sort_values(
            metrics[0], ascending=False).iloc[:n_blend, :]
        self.log(
            'More info about picked checkpoints can be found here: {}/inf_ckpt_selected.csv'.format(self.logs_dir))
        logs_pd_selected.to_csv(self.logs_dir / 'inf_ckpt_selected.csv')
        return logs_pd_selected['id'], logs_pd_selected['config']

    def high_fidelity_training(self, top_n, n_seed):
        self.log('Run configurations for more epochs to achieve better results...')
        logs_summary = auto_summarize_logs(self.dataset)
        gm = GPUManager()
        procs = []
        seeds = [1, 12, 123, 1234, 2, 4, 6, 8]
        for i in range(top_n):
            config = logs_summary.iloc[i, :]['config']
            self.log('Configuration {}: {} ...'.format(i + 1, config))
            config = eval(config)
            config['epochs'] = 2000
            config['note'] = 'more_epochs_run'
            for seed in seeds[:n_seed]:
                config['seed'] = seed
                config['gpu'] = gm.wait_free_gpu(thre=0.5)
                cmd = config2cmd(config)
                p = subprocess.Popen(cmd, shell=True)
                procs.append(p)
                time.sleep(10)
        for p in procs:
            p.wait()
        self.log('Run Complete!', with_time=True)

    def log(self, msg=None, with_time=False):
        msg = str(msg)
        if with_time:
            msg = msg + ' time elapsed {:.2f} hrs ({:.1f} mins)'.format(
                (time.time() - self.start) / 3600.,
                (time.time() - self.start) / 60.
            )
        with open(self.logs_dir / 'inference_log.txt', 'a+') as f:
            f.write(msg + '\n')
            print(msg)

    def pasp(self):
        for level in [1, 2, 3]:
            self.log('Run solution for perturbed test level {}...'.format(level))
            M, M_prime, Q, Q_prime = perturb_test(self.dataset, level)

            P = self.blend_and_inference(custom_dataset=M)
            P_prime = self.blend_and_inference(custom_dataset=M_prime)  # P'

            L_P_P_prime = regression_metrics(P.numpy(), y_pred=P_prime.numpy())
            L_Q_Q_prime = regression_metrics(Q, y_pred=Q_prime)
            self.log('L(P, P\') is {}, and\n L(Q, Q\') is {}'.format(
                L_P_P_prime, L_Q_Q_prime))
            self.log('\Delta_RMSE={}'.format(
                L_P_P_prime['rmse'] - L_Q_Q_prime['rmse']))