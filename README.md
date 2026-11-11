# directed-network-simulations

有向ネットワーク上での SAR（Susceptible-Adopted-Recovered）モデルによる情報拡散・感染症伝播シミュレーションを行う Java プロジェクトです。

## 概要

本プロジェクトは、有向辺と無向辺が混在するネットワーク上で、閾値付き SAR モデルのイベント駆動シミュレーションを高速に実行します。複数のネットワークトポロジモデル（Configuration Model、Barabási-Albert モデル、Erdős-Rényi モデル、実ネットワーク等）をサポートし、ForkJoinPool による並列バッチ処理で大規模なパラメータスキャンを効率的に行えます。

C 言語による EBCM（Edge-Based Compartmental Model）数値計算と、Python/Jupyter による EBCM 解析も含みます。

## 必要環境

### Java シミュレーション

- **Java 21** 以上
- **Gradle 9.x**（Gradle Wrapper 同梱のため別途インストール不要）

### C EBCM

- **Clang**（C11 対応）
- **OpenMP**（`same-in-out.c` の並列実行に使用）

### Python EBCM（Jupyter ノートブック）

- **Python 3.x**
- **Jupyter Notebook / JupyterLab**
- **NumPy**, **SciPy**, **Matplotlib**

---

## セットアップ

### Java 21 のインストール

#### macOS（Homebrew）

```bash
brew install --cask temurin@21
```

#### Linux（SDKMAN）

```bash
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk install java 21-tem
```

インストール後、バージョンを確認します。

```bash
java -version
```

### OpenMP のインストール（macOS）

macOS の Clang は標準で OpenMP を含まないため、`libomp` を別途インストールします。

```bash
brew install libomp
```

### プロジェクトのセットアップ

```bash
git clone <repository-url>
cd directed-network-simulations
./gradlew build
```

---

## 実行方法

### SAR シミュレーション

```bash
./gradlew :app:runSAR
```

`sim.SAR` クラスの `SimulationConfig` を編集することで、ネットワークタイプ・頂点数・感染率・閾値などのパラメータを変更できます。

### グラフ生成

```bash
./gradlew :app:runGraphGen
```

グラフを生成し、NetworkX 互換のエッジリスト形式で `out/edgelist/` に出力します。`SimulationConfig` の `loadFromEdgeList = true` と組み合わせて使います。

### C EBCM

```bash
bash scripts/run-c.sh
```

内部では以下のようなコマンドを実行しています。

```bash
clang -fopenmp -std=c11 -Wall -Wextra -Wpedantic -O0 -g c-lang/ebcm/same-in-out.c -o build/same-in-out
./build/same-in-out
```

### テスト

```bash
./gradlew test
```

---

## プロジェクト構成

```
app/src/main/java/sim/
├── App.java                          # デモ・動作確認用エントリポイント
├── SAR.java                          # SAR シミュレーション（並列バッチ処理）
├── GraphGen.java                     # グラフ生成・エッジリスト出力
├── network/
│   ├── DirectedGraph.java            # 有向グラフ（CSR 形式）
│   └── topology/
│       ├── DirectedBA.java           # Barabási-Albert モデル（有向）
│       ├── DirectedCM.java           # Configuration Model（基本）
│       ├── DirectedCMOutPow.java     # CM（出次数パワーロー分布）
│       ├── DirectedCMInPow.java      # CM（入次数パワーロー分布）
│       ├── PowPow.java               # CM（入出次数ともパワーロー分布）
│       ├── SameInOut.java            # CM（入次数＝出次数、パワーロー分布）
│       ├── EgoTwitter.java           # ego-Twitter 実ネットワーク（SNAP）
│       ├── RevEgoTwitter.java        # ego-Twitter の辺を反転
│       ├── Gplus.java                # Google+ 実ネットワーク（SNAP）
│       ├── RevGplus.java             # Google+ の辺を反転
│       ├── HiggsSocial.java          # Higgs Twitter 実ネットワーク（SNAP）
│       ├── RevHiggsSocial.java       # Higgs Twitter の辺を反転
│       └── undirected/
│           ├── BA.java               # Barabási-Albert モデル（無向）
│           ├── CM.java               # Configuration Model（無向）
│           └── ER.java               # Erdős-Rényi モデル
├── simulation/
│   ├── SARSimulator.java             # イベント駆動 SAR シミュレータ
│   ├── SARGillespieSimulator.java    # Gillespie 法 SAR シミュレータ
│   └── SARResult.java                # シミュレーション結果・CSV 出力
└── utils/
    ├── ArrayUtils.java               # 配列ユーティリティ
    ├── RandomUtils.java              # 確率分布ユーティリティ
    ├── PathsEx.java                  # パスユーティリティ
    └── SwitchUtils.java              # ネットワークタイプのディスパッチ

c-lang/ebcm/                          # C 言語 EBCM 数値計算
├── directed-infty.c                  # 有向ネットワーク（最大次数 ∞）
├── same-in-out.c                     # 入次数＝出次数（OpenMP 並列）
└── undirected-infty.c                # 無向ネットワーク（最大次数 ∞）

notebooks/                            # Jupyter ノートブック（Python EBCM）
├── directed-ebcm.ipynb
└── utils/
    ├── directed_ebcm.py
    ├── ebcm.py
    ├── calculation.py
    └── switch.py
```

---

## DirectedGraph クラスの使い方

`DirectedGraph`（`sim.network.DirectedGraph`）は本プロジェクトの中核となるクラスで、**有向辺と無向辺が混在するグラフ**を CSR（Compressed Sparse Row）形式で効率的に表現します。

### 基本概念

本プロジェクトのグラフでは、各辺が**有向辺（directed）**と**無向辺（nondirected）**のいずれかに分類されます。

- **有向辺** `u → v`: 頂点 `u` から `v` への一方向の辺
- **無向辺** `u -- v`: 内部的には `u → v` と `v → u` の 2 本の有向辺として格納され、`isUndirected` フラグで無向由来であることが識別されます

この区別は SAR シミュレーションにおいて重要です。有向辺と無向辺で異なる感染率（`lambdaDirected`, `lambdaNondirected`）を設定できます。

### 主要プロパティ

| プロパティ | 型 | 説明 |
|---|---|---|
| `name` | `String` | グラフの識別名（例: `"DirectedBA"`, `"DirectedCMOutPow"`） |
| `n` | `int` | 頂点数 |
| `m` | `int` | 展開後の有向辺の総数（無向辺は 2 本としてカウント） |

### グラフの生成

`DirectedGraph` はプライベートコンストラクタを持ち、以下の方法で生成します。

#### 1. トポロジモデルから生成する（推奨）

各トポロジモデルの `generate()` メソッドを使います。

```java
import sim.network.DirectedGraph;
import sim.network.topology.*;
import sim.network.topology.undirected.*;

// --- Barabási-Albert モデル（有向）---
// N=10000, 初期完全グラフ6頂点, 新規ノードあたり5弧, シード42
DirectedGraph baDirected = DirectedBA.generate("DirectedBA", 10000, 6, 5, 42L);

// --- Barabási-Albert モデル（無向）---
DirectedGraph baUndirected = BA.generate("BA", 10000, 6, 5, 42L);

// --- Configuration Model（出次数パワーロー分布） ---
// N=100000, 最小出次数5, 最大出次数316, 平均無向次数0, 指数2.5, シード42
DirectedGraph cmOutPow = DirectedCMOutPow.generate(
    "DirectedCMOutPow", 100_000, 5, 316, 0.0, 2.5, 42L
);

// --- Configuration Model（入次数パワーロー分布） ---
// 内部的に DirectedCMOutPow で生成後、有向辺を反転する
DirectedGraph cmInPow = DirectedCMInPow.generate(
    "DirectedCMInPow", 100_000, 5, 316, 0.0, 2.5, 42L
);

// --- Configuration Model（基本） ---
// N=1000, 平均次数10, シード1234567890
DirectedGraph cmBasic = DirectedCM.generate("DirectedCM", 1000, 10, 1234567890L);

// --- Erdős-Rényi モデル（全辺無向） ---
// N=10000, 平均次数6.0, シード42
DirectedGraph erGraph = ER.generateERFromKAve(10000, 6.0, 42L);

// --- Erdős-Rényi モデル（有向） ---
// 期待平均入次数・出次数6.0、シード42
DirectedGraph directedErGraph = DirectedER.generateFromMeanDegree(10000, 6.0, 42L);
```

#### 2. エッジリストから直接生成する

辺のソース・デスティネーション・無向フラグの配列から生成できます。

```java
int n = 5; // 頂点数
int[] srcs =         {0, 1, 2, 3};
int[] dsts =         {1, 2, 3, 4};
boolean[] isUndirected = {false, true, false, true};
// 辺: 0→1 (有向), 1--2 (無向), 2→3 (有向), 3--4 (無向)

DirectedGraph g = DirectedGraph.fromEdgeListWithUndirectedFlag(
    "myGraph", n, srcs, dsts, isUndirected
);
// 内部的には無向辺が展開され、m = 4 + 2 = 6 本の有向辺として格納される
```

> **注意**: 無向辺 `u -- v` は内部的に `u → v` と `v → u` の 2 本に展開されるため、`m` は元の辺数よりも大きくなります。

### 隣接頂点の走査

CSR 形式のため、特定の頂点の隣接頂点を効率的に列挙できます。

```java
DirectedGraph g = DirectedBA.generate("DirectedBA", 1000, 5, 3, 42L);

int u = 0; // 頂点 0 について調べる

// --- 出隣接（out-neighbors）の走査 ---
DirectedGraph.IntRange outRange = g.outNeighborRange(u);
for (int i = outRange.start; i < outRange.end; i++) {
    int neighbor = g.getOutNeighbor(i);
    boolean isUndirected = g.isOutUndirected(i);
    System.out.println("u -> " + neighbor + (isUndirected ? " (無向)" : " (有向)"));
}

// --- 入隣接（in-neighbors）の走査 ---
DirectedGraph.IntRange inRange = g.inNeighborRange(u);
for (int i = inRange.start; i < inRange.end; i++) {
    int neighbor = g.getInNeighbor(i);
    boolean isUndirected = g.isInUndirected(i);
    System.out.println(neighbor + " -> u" + (isUndirected ? " (無向)" : " (有向)"));
}
```

`IntRange` は `[start, end)` の半開区間を表すクラスで、CSR 形式のポインタ配列に対応しています。

### グラフの変換

```java
// 有向辺の向きを全て逆転（無向辺はそのまま）
DirectedGraph reversed = g.reverseDirectedEdges();

// 名前を指定して逆転
DirectedGraph reversedNamed = g.reverseDirectedEdges("myReversedGraph");
```

`DirectedCMInPow` はこの機能を活用し、`DirectedCMOutPow` で生成したグラフの有向辺を反転することで入次数パワーロー分布を実現しています。

### グラフ情報の表示

```java
g.printInfo();
```

以下のような情報が標準出力に表示されます:

```
=== Graph Info ===
Name: DirectedCMOutPow
Vertices (n): 100000
Directed Arcs (m): 847132
Undirected Arcs (m): 0

Average in-degree (excluding undirected edges): 8.47132
Average in-degree (including undirected edges): 8.47132
Average out-degree (excluding undirected edges): 8.47132
Average out-degree (including undirected edges): 8.47132

Max in-degree: 28
Min in-degree: 0
Max out-degree: 295
Min out-degree: 5

Largest WCC size (directed-only edges): 99996
Largest WCC size (including undirected edges): 99996
```

### 連結性の確認

```java
// 有向辺のみで最大弱連結成分のサイズを取得
int wccDirectedOnly = g.checkConnected(false);

// 有向辺＋無向辺で最大弱連結成分のサイズを取得
int wccAll = g.checkConnected(true);
```

### エッジリストのエクスポート

Python の NetworkX で読み込める形式でエッジリストをファイルに出力できます。

```java
import java.nio.file.Paths;

g.writeEdgeList(Paths.get("output/edgelist.txt"));
```

出力形式（スペース区切り）:
```
0 1
0 3
1 2
...
```

Python での読み込み:
```python
import networkx as nx
G = nx.read_edgelist("output/edgelist.txt", nodetype=int, create_using=nx.DiGraph)
```

---

## SAR シミュレーション

SAR モデルは、ネットワーク上での情報拡散・感染症伝播をシミュレートする閾値付きモデルです。

### ノードの状態遷移

```
S (Susceptible) → A (Adopted) → R (Recovered)
   感受性状態        採用状態        回復状態
```

- **S → A**: ノードが累計で閾値回数の感染を受けると採用状態に遷移
- **A → R**: 指数分布に従う時間経過後に回復状態に遷移

### シミュレータの選択

`SARSimulator`（イベント駆動方式）と `SARGillespieSimulator`（Gillespie 法）の 2 種類が用意されています。`SimulationConfig.useGillespie` フラグで切り替えます。両者は統計的に同値なモデルです。

### 基本的な使い方

```java
import sim.network.DirectedGraph;
import sim.network.topology.DirectedCMOutPow;
import sim.simulation.SARSimulator;
import sim.simulation.SARResult;
import java.util.Arrays;

// 1. グラフを生成
DirectedGraph g = DirectedCMOutPow.generate(
    "DirectedCMOutPow", 100_000, 5, 316, 0.0, 2.5, 42L
);

// 2. パラメータを設定
double lambdaDirected = 0.5;     // 有向辺の感染率
double lambdaNondirected = 0.0;  // 無向辺の感染率
double mu = 1.0;                 // 回復率
double tMax = 200.0;             // シミュレーション最大時刻
int[] thresholdList = new int[g.n];
Arrays.fill(thresholdList, 4);   // 全ノードの閾値を 4 に設定

// 3. 初期感染者をランダムに選択（全体の 10%）
int initialCount = (int) (g.n * 0.1);
int[] initialInfecteds = new int[initialCount];
for (int i = 0; i < initialCount; i++) {
    initialInfecteds[i] = i;
}

// 4. シミュレーション実行
long seed = 12345L;
SARResult result = SARSimulator.simulate(
    g, lambdaDirected, lambdaNondirected, mu, tMax,
    thresholdList, initialInfecteds, seed
);

// 5. 結果の確認
int finalA = result.A.get(result.A.size() - 1);
int finalR = result.R.get(result.R.size() - 1);
System.out.println("最終採用数: " + finalA);
System.out.println("最終回復数: " + finalR);
```

### 結果の CSV 出力

```java
import java.nio.file.Paths;

// 時系列データを出力（各イベント発生時の S, A, R, Phi）
result.writeTimeSeriesCsv(
    Paths.get("out/timeseries.csv"),
    0,                   // イテレーション番号
    0.1,                 // 初期感染率 rho0
    lambdaDirected,
    lambdaNondirected,
    mu,
    false                // 追記モード
);

// 最終状態のみ出力
result.writeFinalStateCsv(
    Paths.get("out/final.csv"),
    0, 0.1, lambdaDirected, lambdaNondirected, mu,
    false
);
```

**時系列 CSV カラム**: `itr, rho_0, lambda_d, lambda_u, mu, time, A, R, Phi`

**最終状態 CSV カラム**: `itr, rho_0, lambda_d, lambda_u, mu, time, initial_adopted_time, final_adopted_time, A, R, Phi`

---

## ネットワークトポロジモデル

### Barabási-Albert モデル（`DirectedBA` / `BA`）

優先的選択（preferential attachment）に基づくスケールフリーネットワークを生成します。

- **有向（`DirectedBA`）**: 新規ノードは既存ノードの**入次数**に比例する確率で選択
- **無向（`BA`）**: 新規ノードは既存ノードの**次数**に比例する確率で選択

```java
DirectedGraph gDirected   = DirectedBA.generate("DirectedBA", 10000, 6, 5, 42L);
DirectedGraph gUndirected = BA.generate("BA", 10000, 6, 5, 42L);
```

### Configuration Model - 出次数パワーロー（`DirectedCMOutPow`）

出次数がパワーロー分布 `P(k) ~ k^{-gamma}` に従い、入次数は総数をランダムに配分、無向次数はポアソン分布に従うグラフを生成します。

```java
DirectedGraph g = DirectedCMOutPow.generate(
    "DirectedCMOutPow",
    100_000,  // 頂点数
    5,        // 最小出次数 kOutMin
    316,      // 最大出次数 kOutMax
    0.0,      // 平均無向次数 kuAve
    2.5,      // パワーロー指数 gamma
    42L       // シード
);
```

### Configuration Model - 入次数パワーロー（`DirectedCMInPow`）

`DirectedCMOutPow` で生成したグラフの有向辺を反転し、入次数がパワーロー分布に従うグラフを得ます。

```java
DirectedGraph g = DirectedCMInPow.generate(
    "DirectedCMInPow", 100_000, 5, 316, 0.0, 2.5, 42L
);
```

### Configuration Model - 入出次数ともパワーロー（`PowPow`）

入次数・出次数ともにパワーロー分布に従うグラフを生成します。

```java
DirectedGraph g = PowPow.generate("PowPow", 100_000, 5, 1000, 2.5, 0, 42L);
```

### Configuration Model - 入次数＝出次数（`SameInOut`）

全ノードで入次数と出次数が等しく、パワーロー分布に従うグラフを生成します。

```java
DirectedGraph g = SameInOut.generate("SameInOut", 100_000, 5, 1000, 2.5, 42L);
```

### Configuration Model - 基本（`DirectedCM`）

平均次数 `kHat` を指定して、有向辺と無向辺が混在するグラフを生成します。入次数は全ノード一律 `kHat`、出次数は `[0, 2*kHat]` の一様分布で総和が `n*kHat` になるよう調整、無向次数は `2*kHat - ko` で決定されます。

```java
DirectedGraph g = DirectedCM.generate("DirectedCM", 1000, 10, 42L);
```

### Erdős-Rényi モデル（`ER`）

古典的なランダムグラフモデルです。**全ての辺が無向**として生成されます。

```java
// 確率指定
DirectedGraph gFromP = ER.generateERFromP(10000, 0.001, 42L);

// 平均次数指定
DirectedGraph gFromK = ER.generateERFromKAve(10000, 6.0, 42L);
```

### 有向 Erdős-Rényi モデル（`DirectedER`）

自己ループを除く各順序付き頂点対 `u -> v` を独立に確率 `p` で結びます。
自己ループと多重辺は生成しません。入次数・出次数は厳密には
`Binomial(N - 1, p)` に従い、疎な場合に平均 `(N - 1) * p` のポアソン分布へ近似します。

```java
// 順序付き頂点対ごとの確率を指定
DirectedGraph gFromP = DirectedER.generateFromP(10000, 0.001, 42L);

// 期待平均入次数および出次数を指定
DirectedGraph gFromK = DirectedER.generateFromMeanDegree(10000, 6.0, 42L);
```

### 実ネットワーク（ego-Twitter / gplus / higgs-social）

[SNAP（Stanford Network Analysis Project）](https://snap.stanford.edu/data/) からダウンロードしたエッジリストを読み込みます。

| ネットワーク名 | ファイルパス |
|---|---|
| ego-Twitter | `app/data/real-network/twitter_combined.txt` |
| gplus | `app/data/real-network/gplus_combined.txt` |
| higgs-social | `app/data/real-network/higgs-social_network.edgelist` |

```java
DirectedGraph twitter = EgoTwitter.loadFromDefaultEdgeList();
```

`Rev*` クラス（`RevEgoTwitter`, `RevGplus`, `RevHiggsSocial`）は対応するグラフの有向辺を反転します。

---

## 出力ディレクトリ構造

SAR シミュレーションの結果は以下の構造で保存されます:

```
out/sar/{optionPath}/threshold={T}/{networkType}/N={N}/{networkSpecific}/results_{batch_idx}.csv
```

例:
```
out/sar/in-out/threshold=3/PowPow/N=500000/gamma=2.50/kdMin=5/kdMax=1000/swapNum=0/results_00.csv
```

## Java内部の責務分担

実験条件は従来どおり各エントリポイントの設定欄で編集します。SAR・グラフ生成・パーコレーションの実行方法と出力形式は維持しています。

- `SwitchUtils.GraphGeneratorParams.builder(type, n)` は名前付きで生成条件を組み立てます。同じ設定を `generateGraph` と `buildNetworkPath` に渡せます。実頂点数でパスを作る場合は `withN` を使います。従来の引数形式も利用できます。
- `SARGraphPreparation` はSARのグラフ準備とファイル探索、`InitialAdopters` は固定シードでの初期感染者選択を担当します。計算方式と記録間隔の選択は `SAR` 内で行います。
- `DirectedGraph` はCSR表現と公開APIを保持し、内部の `GraphOperations` と `FflRewiringState` に指標計算・変換を委譲します。実ネットワークの読み込みは共通のパーサーを使用します。
- `SARResult` は結果データと公開APIを保持し、CSV出力を `SARResultCsvWriter` に委譲します。

固定シードの互換性テストは、整理前のコードから取得した辺の順序と両シミュレーターの時系列・感染時刻のSHA-256を基準にしています。BA設定ヘルパーの引数割り当て修正は、専用の回帰テストで検証します。

## SPIR / SIR シミュレーション

```bash
bash scripts/run-spir.sh
# または
./gradlew :app:runSPIR
```

`sim.SPIR.SimulationConfig` で設定します。`mode = ModelMode.SPIR` では
P（症状なし・感染力あり）を加え、S→P→I→R と進みます。
P→S 辺は `lambda1`、I→S 辺は `lambda2` のレートでSを感染させます。
P→I は `mu1`、I→R は `mu2` です。既定では連続時間のGillespie法で処理し、
有向辺は向きに沿い、無向辺は両方向に同じレートを適用します。

`mode = ModelMode.SIR` では初期感染者も新規感染者も直接Iになり、Pを省略します。
`lambda2` と `mu2` が通常のSIRの感染・回復レートです。
SPIRでは4レートのリストの全組み合わせ、SIRでは `lambda2List × mu2List` だけを実行します。
初期感染人数はSARと同じ `max(1, floor(N × rho0))` で、SPIRは全員P、SIRは全員Iです。
同じバッチ・反復では、レート条件間でオブサーバー・初期ノード・シミュレーションシードを共用します。

ネットワークは `GraphGeneratorParams.builder(networkType, N)` で指定します。
既存SARと同じ生成器・エッジリスト探索・ランダマイズを使用します。
実行前に `SimulationConfig` のネットワーク・レート・バッチ数・反復数を確認してください。
実験条件に合わせて変更してください。

出力先は `app/out/spir/{optionPath}/model={spir|sir}/{networkPath}/randomization={mode}/phi={phi}/`。
Gradleタスクの作業ディレクトリが `app/` のため、実行時はこの場所になります。
バッチ実行のCSV列は `model,itr,rho_0,lambda1,lambda2,mu1,mu2,time,S,P,I,R,phi,observer_count,Q`。
`Q` は累積隔離人数で、健康状態とは別に数えるため `S+P+I+R=N` にQを加算しません。
SIRの `lambda1`・`mu1` 列は空欄です。`isFinal=true` は `tMax` 時点の1行、
`false` は `dt` 間隔の時系列を出力します。終息後の状態も最後の格子点まで記録します。
レート0を許可し、負値・非有限値は拒否します。
実行サマリーJSONには設定・シード・出力先・完了数・エラーを保存します。
既存CSVは上書きせず、番号付きの別ファイルを作成します。

### オブサーバーと隔離

`SPIR.java` 最下部の `SimulationConfig.phi`（有限の0〜1）で割合を指定します。
`phi = 0.0` でオブサーバーなしのモデルになります。
全頂点から `floor(N × phi)` 人を重複なしで無作為抽出し、初期感染者は残りの頂点から抽出します。
初期感染人数は全頂点数に対する `rho0` で決め、オブサーバーを除いた人数を分母にはしません。
非オブサーバーが不足する場合はエラーにします。例えば初期感染者100人を要求する設定では、
全員がオブサーバーになる `phi=1` は実行できません。
抽出には感染シミュレーションと独立の固定シードを使用します。

オブサーバーも通常のS/P/I/R状態を取り、Pの間は感染源になります。
監視対象は隣接頂点で、有向グラフでは入辺・出辺の両方の隣接頂点です。
Pへの感染は検知せず、Iへの遷移を待ち時間なしで検知します。
検知したオブサーバーは、自身と自身の隣接頂点を隔離しますが、すでにIの隣接頂点は除外します。
オブサーバー自身がIになった場合は自身も隔離します。
同じ発症イベントを検知した複数のオブサーバーは、そのイベント時点で未隔離であれば全員が反応します。

隔離は解除せず、監視も含めた接触辺を切断し、感染の送受信を完全に止めます。
隔離されたオブサーバーは監視を停止し、隔離済みPが後でIになっても新たな検知は起こしません。
隔離しても健康状態は変えず、P→I→Rは元の `mu1`・`mu2` で進行します。
このため、隔離済みPの発症後やオブサーバーの自己発症時には、隔離中のIが存在します。
SIRモードでは時刻0の初期Iにも検知を適用します。
`SPIRResult` はオブサーバーフラグ、隔離フラグ、各頂点の隔離時刻も保持します。
既存のオブサーバー引数なしのシミュレーションAPIは、オブサーバーなしとして動作します。

### SPIRの方式選択と性能

`SimulationConfig.engine` で `SPIRSimulator.Engine.GILLESPIE` または
`SPIRSimulator.Engine.EVENT_DRIVEN` を選べます。既定はGillespieです。
両方式は状態遷移・監視・隔離を共有し、各方式内ではシードで再現できます。
方式間や旧実装との軌跡の一致は保証しません。
最終状態のみの実行は途中履歴を保存せず、グラフの読み取り専用索引をバッチ内で再利用します。
既存の `SPIRGillespieSimulator.simulate` は、従来どおりイベント時刻または固定時刻の履歴を返します。
CSV形式と出力先は維持しています。

`SimulationConfig.equalLambda = true` にすると、`lambda2List` の各値を
`lambda1 = lambda2` として走査します。この場合、`lambda1List` は走査に使用しません。
例えば `lambda2List = ArrayUtils.arange(0.0, 4.0, 0.2)` なら、
`(0,0), (0.2,0.2), …, (4,4)` の21組です。
`mu1List × mu2List` の走査は維持します。`false` は従来の全組み合わせ走査です。
SIRモードはこの設定によらず `lambda2List × mu2List` を走査します。

### TwoPointDirectedCM：2点次数分布と同一頂点の入出次数相関

入次数・出次数はともに確率 `p` で `k1`、確率 `1-p` で `k2` を取ります。
`r` は同じ頂点の入次数・出次数の Pearson 相関係数です（隣接頂点間の相関とは異なります）。
`1 <= k1 < k2`、`0 < p < 1`、`-min(p/(1-p), (1-p)/p) <= r <= 1` が必要です。
自己ループ・多重辺を許容し、連結性は保証しません。

```java
var params = SwitchUtils.GraphGeneratorParams.builder("TwoPointDirectedCM", 1000)
        .k1(2).k2(8).p(0.5).r(0.5).build();
var graph = SwitchUtils.generateGraph(params, 42L);
var counts = sim.network.topology.TwoPointDirectedCM.degreePairCounts(1000, 0.5, 0.5);
System.out.println(counts.realizedP() + " " + counts.realizedR());
```

理論同時分布は `a=p*p+r*p*(1-p)` として、
`P11=a, P12=P21=p-a, P22=1-2*p+a` です。
有限サイズでは `n1=floor(N*p+0.5), n2=N-n1, pN=n1/N`、
`b=clip(floor(N*pN*(1-pN)*(1-r)+0.5), 0, min(n1,n2))` とし、
人数を `(n1-b,b,b,n2-b)` に固定します。
両次数クラスの人数が正になるNが必要です。実現相関は `1-N*b/(n1*n2)` です。
seedは頂点への割り当てと接続を変えますが、人数・平均次数は変えません。

SARでは `SimulationConfig` の `networkType` と `k1,k2,p,r` を設定します。
GraphGenにも同じ設定項目があります。相関比較には `randomizationMode=NONE` を使用してください。
入次数列と出次数列を独立にシャッフルする変換は同一頂点の相関を変えます。
出力パスには指定した `k1,k2,p,r` が入り、実現値とは区別します。

ローカル解析用 `test-notebooks/utils/directed_ebcm.py` では次のように構築できます。

```python
from utils.directed_ebcm import build_two_point_joint_degree_dist, build_binom, Q_dir

theory = build_two_point_joint_degree_dist(2, 8, 0.5, 0.5)
finite = build_two_point_joint_degree_dist(2, 8, 0.5, 0.5, n=1000)
# Javaとの比較にはfiniteを使う。既存のDynamicsConfigを指定する。
# Q_dir(finite, dynamics, build_binom(finite.ki_max), eps=0.2)
```

同時分布の行が入次数、列が出次数で、次数値そのものが添字です。
Java出力は既存の `build_joint_degree_dist_from_edge_list` でも読み込めます。
この例は設定方法を示すもので、実験実行は別途条件を確定して行います。
