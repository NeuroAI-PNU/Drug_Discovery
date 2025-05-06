# Drug_Discovery Protocol
## 폴더 설명
### Workflow : Drug discovery를 위한 일련의 과정을 코드와 함께 순서대로 정리한 폴더
### package : 사용되는 코드 중 주요 클래스 및 함수를 모듈화하여 정리한 폴더

## Workflow
#### 0 단계 : 데이터 전처리 단계, ZINC20으로부터 화합물 데이터를 핸들링하고 다음 과정에서 사용될 단백질 선택 단계
 - 0-1 : 다운로드 받은 ZINC20 데이터셋에서 zinc_id 와 SMILES만 파싱
 
   <https://zinc20.docking.org/tranches/home/#> 에서 In-Stock 설정으로 다운로드 받음
   
 -  0-2 : 분자 도킹 시뮬레이션을 통해 타겟 단백질의 구조-참조물질(reference chemical) 사이 구조 친화도를 분석하고 그 결과를 기반으로 하나의 단백질 구조 선택
 -  0-3 : 0-2를 통해 얻은 도킹 점수 정리하는 코드

#### 1 단계 : Lipinski's rule of Five를 통해 Drug-likeness한 화합물만을 선정하는 단계
 - 1-1 : Lipinski 법칙 중 분자량, LogP, Hydrogen Donor, Hydrogen Acceptor 4가지 기준에 부합하는 화합물만 선택
 - 1-2 : Lipinski 분석 결과 정리하는 코드

#### 2 단계 : Pan-assay interference compounds(PAINS) 잔기 제거 단계
 - 2 : PAINS 잔기를 포함하는 화합물을 식별하고 그 화합물은 리스트에서 제거

#### 3 단계 : 참조물질과 ZINC20으로부터 얻은 화합물 사이 구조적 유사도 분석 단계
 - 3-1 : tanimoto similiarty를 통해 두 물질 사이 유사도 점수 계산
 - 3-2 : threshold = 0.8을 기준으로 그 이상의 유사도 점수를 가지는 화합물만 선택
 - 3-3 : 전체 화합물과 0.8 이상 화합물이 가지는 분포도를 보여주는 그래프 형성

#### 4 단계 : 참조물질과 ZINC20으로부터 얻은 화합물 사이 분자특징적 유사도 분석 단계
 - 4-1 : 1-1에서 사용한 4가지 기준을 기반으로 k-means 알고리즘을 이용한 참조물질과 동일한 cluster에 존재하는 화합물 선택
 - 4-2 : 참조물질 그리고 참조물질과 같은 clsuter에 존재하는 화합물 비교 (필요시 이 단계 진행)

#### 5 단계 : tanimoto(3단계)를 통해 얻은 화합물과 k-means(4단계)를 통해 얻은 화합물들 사이 공통 화합물 식별 단계
 - 5 : 3,4단계를 통해 얻은 데이터 정리 및 공통 화합물 식별

#### 6 단계 : 분자 도킹 시뮬레이션 수행 단계
 - 6-1 : tanimoto similairt를 통해 얻은 화합물과 0단계에서 선택된 단백질 구조 사이 구조적 친화도 평가

   이 과정에서 memory issue로 인해 multiprocessing 방법으로 수행함.
 - 6-2 : k-means를 통해 얻은 화합물과 0단계에서 선택된 단백질 구조 사이 구조적 친화도 평가

#### 7 단계 : 도킹 점수 parsing 단계
 - 7-1 : 6-1를 통해 얻은 도킹 점수 parsing
 - 7-2 : 6-2를 통해 얻은 도킹 점수 parsing

#### 8 단계 : parsing된 도킹 점수 분석 단계
 - 8-1 : 참조물질의 도킹 점수와 함께 7 단계에서 얻은 도킹 점수의 분포도를 하나의 그래프에 나타내는 코드

#### 9 단계 : AI(SOTA 모델 + RDKit)를 이용한 물리화학적 특징 예측 단계
 - 9-1 (폴더) : Blood-Brain Barrier Permeability 분석 (model : GLAM)
 - 9-2 (폴더) : Biological Toxicity 분석 (model : elEmBERT)
 - 9-3 : 합성용이도를 평가 (model : SAScoerer)

   SOTA model reference
   <https://github.com/yvquanli/GLAM>

   <https://github.com/dmamur/elembert>

   <https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py>
