# Drug_Discovery Protocol
## 폴더 설명
### Workflow : 구조&화학적특징 유사도 + AI 기반 Virtual Screening을 위한 코드 및 일련의 순서를 정리해둔 폴더
### MolecularDynamics : Virtual Screening 후 검증을 위한 분자동역학 시뮬레이션 코드를 정리해둔 폴더
### package : 사용되는 코드 중 주요 클래스 및 함수를 모듈화하여 정리한 폴더


## Workflow
#### 0 단계 : 데이터 전처리 단계, ZINC20으로부터 화합물 데이터를 핸들링하고 다음 과정에서 사용될 단백질 선택 단계
 - 0-1 : 다운로드 받은 ZINC20 데이터셋에서 zinc_id 와 SMILES만 파싱
 
   <https://zinc20.docking.org/tranches/home/#> 에서 In-Stock 설정으로 다운로드 받음
   
 -  0-2 : 분자 도킹 시뮬레이션을 통해 타겟 단백질의 구조-참조물질(reference chemical) 사이 구조 친화도를 분석하고 그 결과를 기반으로 하나의 단백질 구조 선택
 -  0-3 : 0-2를 통해 얻은 도킹 점수 정리하는 코드

#### 1 단계 : Lipinski's rule of Five를 통해 Drug-likeness한 화합물만을 선정하는 단계
 - 1-1 : Lipinski 법칙 중 분자량, LogP, Hydrogen Donor, Hydrogen Acceptor 4가지 기준에 부합하는 화합물만 선택.
 - 1-2 : Lipinski 분석 결과 정리하는 코드.

#### 2 단계 : Pan-assay interference compounds(PAINS) 잔기 제거 단계
 - 2 : PAINS 잔기를 포함하는 화합물을 식별하고 그 화합물은 리스트에서 제거.

#### 3 단계 : 참조물질과 ZINC20으로부터 얻은 화합물 사이 구조적 유사도 분석 단계
 - 3-1 : tanimoto similiarty를 통해 두 물질 사이 유사도 점수 계산.
 - 3-2 : threshold = 0.8을 기준으로 그 이상의 유사도 점수를 가지는 화합물만 선택.
 - 3-3 : 전체 화합물과 0.8 이상 화합물이 가지는 분포도를 보여주는 그래프 형성.

#### 4 단계 : 참조물질과 ZINC20으로부터 얻은 화합물 사이 분자특징적 유사도 분석 단계
 - 4-1 : 1-1에서 사용한 4가지 기준을 기반으로 k-means 알고리즘을 이용한 참조물질과 동일한 cluster에 존재하는 화합물 선택.
 - 4-2 : 참조물질 그리고 참조물질과 같은 clsuter에 존재하는 화합물 비교. (필요시 이 단계 진행)

#### 5 단계 : tanimoto(3단계)를 통해 얻은 화합물과 k-means(4단계)를 통해 얻은 화합물들 사이 공통 화합물 식별 단계
 - 5 : 3,4단계를 통해 얻은 데이터 정리 및 공통 화합물 식별.

#### 6 단계 : 분자 도킹 시뮬레이션 수행 단계
 - 6-1 : tanimoto similairt를 통해 얻은 화합물과 0단계에서 선택된 단백질 구조 사이 구조적 친화도 평가.

   이 과정에서 memory issue로 인해 multiprocessing 방법으로 수행함.
 - 6-2 : k-means를 통해 얻은 화합물과 0단계에서 선택된 단백질 구조 사이 구조적 친화도 평가.

#### 7 단계 : 도킹 점수 parsing 단계
 - 7-1 : 6-1를 통해 얻은 도킹 점수 parsing.
 - 7-2 : 6-2를 통해 얻은 도킹 점수 parsing.

#### 8 단계 : parsing된 도킹 점수 분석 단계
 - 8-1 : 참조물질의 도킹 점수와 함께 7 단계에서 얻은 도킹 점수의 분포도를 하나의 그래프에 나타내는 코드.

#### 9 단계 : AI(SOTA 모델 + RDKit)를 이용한 물리화학적 특징 예측 단계
 - 9-1 (폴더) : Blood-Brain Barrier Permeability 분석. (model : GLAM)
 - 9-2 (폴더) : Biological Toxicity 분석. (model : elEmBERT)
 - 9-3 : 합성용이도를 평가. (model : SAScoerer)

   SOTA model reference
   
   <https://github.com/yvquanli/GLAM>

   <https://github.com/dmamur/elembert>

   <https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py>

#### 10 단계 : AI 예측 결과 통합 단계
 - 10 : 3가지 서로 다른 점수를 하나의 DataFrame으로 합침.

   BBBP -> 0과 1사이 값이며 1에 가까울수록 BBB 투과율이 높게 평가됨.
   
   Toxicity -> 0과 1사이 값이며 1에 가까울수록 독성이 없다고 평가됨. 단, 결과값이 총 12개 항목으로 나타남.
   
   SAScore -> 1과 12사이 값이며 1에 가까울수록 합성하기 용이하다고 평가됨.

   BBBP 결과는 그대로 유지, Toxicity 결과는 각 항목의 점수를 12로 나누고 모두 더하여 하나의 독성 점수로 사용, SAScore는 MinMaxScaler를 이용하여 0과 1사이 값으로 나타냄.

   BBBP와 Toxicity 12가지 항목 중 0.5 기준으로 0.5 미만인 화합물은 모두 제거. SAscore는 0.5 기준으로 0.5 초과인 화합물은 모두 제거.

   Final Score = BBBP + Toxicity - SAScore 식을 통해 점수를 산출하고 오름차순으로 sorting하여 상위 10개 후보물질 식별

#### 11 단계 : AI 예측 결과 그래프 분석 단계
 - 11 : 세 결과를 통합하여 0.5를 기준으로 세가지 모두 충족하는 화합물과 충족하지 못하는 화합물을 하나의 그래프에 나타냄

#### 12 단계 : 상위 후보물질 10개와 참조물질를 비교하기 위한 데이터 통합 단계
 - 12 : 상위 후보물질 10개가 초기 어떤 그룹에서부터 선별된 화합물인지 확인(Tanimoto인지 또는 K-means에 의해 선별 되었는지 확인)하고 참조물질을 DataFrame 마지막에 위치시켜 비교하도록 데이터 구축.

0 ~ 12 단계 : 
-----
