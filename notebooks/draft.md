# 主成分分析論文まとめ

## タイトル

Jolliffe IT, Cadima J. 2016 *Principal component analysis: a review and recent developments.* Phil. Trans. R. Soc. A 374: 20150202.
<br/>
<br/>

## 1. Introduction

### PCA(Principal Component Analysis)

* 多次元データが持つ情報をで出来るだけ保存して低次元空間に情報を縮約するアルゴリズム
* 多次元データを2次元・3次元などに押し込めればデータを**可視化**することが出来る
* 事前の分布の仮定などはなくデータ依存, 説明的(descriptive)な方法
<br/>
<br/>

### 情報量と分散

* 射影したデータのばらつき（variability）が大きいほど元のデータの情報が保存されていると考えられる
<br/>
<br/>
<img src="../data/pca_figure1.jpeg" alt="PCA" width=500 />
<br/>
<br/>

### 目的
* 元のデータを射影したときに分散が最大となるような軸を探したい
<br/>;
<br/>
<br/>

## 2. The basic method

### PCAのキソ
$\thinspace p \thinspace$個の$\thinspace n \thinspace$次元のベクトル $\boldsymbol{x}_1, ...,\boldsymbol{x}_p$ を並べたデータ行列 $\boldsymbol{X}$ を考え, $\boldsymbol{X}$ の各列の線形結合

<br/>

$$ \boldsymbol{X}\boldsymbol{a} = \sum^{p}_{j=1}a_{j}\boldsymbol{x}_j$$

<br/>

の中で最大の分散を与えるものを探す. 分散は以下の式で与えられる.

<br/>

$$var( \boldsymbol{X}\boldsymbol{a}) = \boldsymbol{a'Sa} = (\boldsymbol{a'X'Xa})$$

<br/>

知りたいのは軸の方向だけなので $||\boldsymbol{a}||^2 = 1$ の制約を課す. ラグランジュ関数

<br/>

$$L(\boldsymbol{a}, \lambda) = \boldsymbol{a'Sa} - \lambda(\boldsymbol{a'a}-1)$$

<br/>

を $\boldsymbol{a}$　に関して微分すると

<br/>

$$\boldsymbol{Sa} - \lambda\boldsymbol{a} = \boldsymbol{0} \Leftrightarrow \boldsymbol{Sa} = \lambda\boldsymbol{a}$$

<br/>

が得られる. すなわち $\boldsymbol{a}$ は $\boldsymbol{S}$ の固有ベクトルであり, 分散最大化問題が固有値問題に帰着された.

---
#### 復習：線形代数

任意の $p \times p$ 実対称行列（例えば上記の $\boldsymbol{S}$） は全部で $p$ 個の実数の固有値を持ち, これらのうち異なる固有値に属する固有ベクトルは互いに直交する. また, 実対称行列は適当な直交行列 $\boldsymbol{P}$ を用いて対角化することが出来る.

---
$\boldsymbol{S}$ の異なる固有値に属する固有ベクトルは互いに直交することと, ノルムの制約から, $\boldsymbol{a_{k}}, \boldsymbol{a_{k'}}$ をそれぞれ $\boldsymbol{S}$ の固有ベクトルとするとき, 

<br/>

```math
\boldsymbol{a_{k}'a_{k'}} = \left\{
\begin{array}{ll}
1 & (k = k') \\
0 & (otherwise)
\end{array}
\right.
```

<br/>

である. これより, 各固有ベクトルによって得られる $\boldsymbol{X}$ の列ベクトルの線形結合　

<br/>

$$\boldsymbol{X}\boldsymbol{a}_k = \sum^{p}_{j=1}a_{jk}\boldsymbol{x}_j$$

<br/>

において, $k \neq k'$ のとき, これらの内積は

<br/>

$$\boldsymbol{a}_{k'}'\boldsymbol{X}' \boldsymbol{X}\boldsymbol{a}_k = \lambda_{k}\boldsymbol{a_{k}'a_{k'}} = 0$$

<br/>

より, 無相関になっていることが分かる.　下図がイメージ. 

<br/>
<br/>
<img src="../data/pca_2.png" alt="PCA_2" width=500 />
(https://statistics.co.jp/reference/software_R/statR_9_principal.pdf   より引用)
<br/>
<br/>

この時の各線形結合 $\boldsymbol{X}\boldsymbol{a}_k$ を**主成分**(Principal Component), 固有ベクトル $\boldsymbol{a}_k$ の各成分を**主成分負荷量**, $\boldsymbol{X}\boldsymbol{a}_k$ の各成分を**主成分得点**と呼ぶ.

<br/>

### 特異値分解（Singular Value Decomposition）

次元の削減を行う手法という意味では主成分分析と似たような手法. 

任意の ランクが $r$ である $n \times p$ 行列 $\boldsymbol{Y}$ は 

<br/>

$$\boldsymbol{Y} = \boldsymbol{ULA'}$$

<br/>

ここで, $\boldsymbol{U}$, $\boldsymbol{A}$ はそれぞれ $n \times r, p \times r$ 行列で各列は直交している. $\boldsymbol{L}$ は $r \times r$ の対称行列であり, その対角成分は $\boldsymbol{Y}$ の**特異値**と呼ばれる. $\boldsymbol{L} = diag(\sigma_1, ..., \sigma_r), \thinspace \sigma_1 \geq \sigma_2 \geq ... \geq \sigma_r \geq 0$ とするとこの分解は一意に定まる. 下図が分解のイメージ. (文字が違うので雰囲気だけ)

<br/>
<br/>
<img src="../data/svd.png" alt="svd" />
(https://qiita.com/kidaufo/items/0f3da4ca4e19dc0e987e   より引用)
<br/>
<br/>

### 特異値分解の性質

上の特異値分解の結果を使い, $q < r$ に対して $\boldsymbol{U}$ の列のうち左から $q$ 本を抜き出した行列を $\boldsymbol{U}_q$, $\boldsymbol{A}$ の列のうち左から $q$ 本を抜き出した行列を $\boldsymbol{A}_q$, $\sigma_1, ..., \sigma_q$ を対角成分に持つ行列を $\boldsymbol{L}_q$ として、

<br/>

$$\boldsymbol{Y}_q = \boldsymbol{U}_q\boldsymbol{L}_q\boldsymbol{A}_q'$$

<br/>

とすると, $\boldsymbol{Y}_q$ はランク $q$ の行列のうち, $\boldsymbol{Y} - \boldsymbol{Y}_q$ の各成分の2乗和を最小にする行列, すなわち成分の2乗誤差の意味で　$\boldsymbol{Y}$ を最もよく近似する行列になっている. 行列 $\boldsymbol{Y}$ を良く近似するより低いランクの行列を求めることは**低ランク近似**と呼ばれる.

<br/>

### 近似の評価

ランク $q$ の近似がどれだけ良いものかは以下の式で測ることが出来る

<br/>

$$\pi = \frac{\sum_{i=1}^{q}\lambda_i}{\sum_{j=1}^{p}\lambda_j}$$

<br/>

この値が70％を超えるような $q$ を選ぶのがcommon practiceであると論文中では述べられている.

<br/>

### 実践：特異値分解

行列の低ランク近似を実際にやってみる. 近似の様子がわかりやすいように,  <a href="http://cbcl.mit.edu/software-datasets/heisele/facerecognition-database.html">The MIT-CBCL face recognition database</a> からダウンロードしてきた顔画像を用いる.

<br/>
<br/>
<div align="center">
    <img src="../data/face.png" alt="svd" width=350/>
    <div>(使用する顔画像)</div>
</div>
<br/>
<br/>
