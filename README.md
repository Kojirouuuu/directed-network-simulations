# directed-network-simulations

有向ネットワーク上での SAR（Susceptible-Adopted-Recovered）モデルによる情報拡散・感染症伝播シミュレーションを行う Java プロジェクトです。

## 概要

本プロジェクトは、有向辺と無向辺が混在するネットワーク上で、閾値付き SAR モデルのイベント駆動シミュレーションを高速に実行します。複数のネットワークトポロジモデル（Configuration Model、Barabási-Albert モデル、Erdős-Rényi モデル）をサポートし、ForkJoinPool による並列バッチ処理で大規模なパラメータスキャンを効率的に行えます。

## 必要環境

- **Java 21** 以上
- **Gradle 9.x**（Gradle Wrapper 同梱）

## セットアップ

```bash
git clone <repository-url>
cd directed-network-simulations
./gradlew build
```

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

BA モデルのグラフを生成し、NetworkX 互換のエッジリスト形式で出力します。

### テスト

```bash
./gradlew test
```

## プロジェクト構成

```
app/src/main/java/sim/
├── App.java                          # デモ・動作確認用エントリポイント
├── SAR.java                          # SAR シミュレーション（並列バッチ処理）
├── GraphGen.java                     # グラフ生成ユーティリティ
├── network/
│   ├── DirectedGraph.java            # 有向グラフ（CSR 形式）
│   └── topology/
│       ├── BA.java                   # Barabási-Albert モデル
│       ├── DirectedCM.java           # Configuration Model（基本）
│       ├── DirectedCMOutPow.java     # CM（出次数パワーロー分布）
│       ├── DirectedCMInPow.java      # CM（入次数パワーロー分布）
│       └── undirected/
│           └── ER.java               # Erdős-Rényi モデル
├── simulation/
│   ├── SARSimulator.java             # イベント駆動 SAR シミュレータ
│   └── SARResult.java                # シミュレーション結果・CSV 出力
└── utils/
    ├── ArrayUtils.java               # 配列ユーティリティ
    ├── RandomUtils.java              # 確率分布ユーティリティ
    └── PathsEx.java                  # パスユーティリティ
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
| `name` | `String` | グラフの識別名（例: `"BA"`, `"DirectedCMOutPow"`） |
| `n` | `int` | 頂点数 |
| `m` | `int` | 展開後の有向辺の総数（無向辺は 2 本としてカウント） |

### グラフの生成

`DirectedGraph` はプライベートコンストラクタを持ち、以下の方法で生成します。

#### 1. トポロジモデルから生成する（推奨）

各トポロジモデルの `generate()` メソッドを使います。

```java
import sim.network.DirectedGraph;
import sim.network.topology.*;
import sim.network.topology.undirected.ER;

// --- Barabási-Albert モデル ---
// N=10000, 初期完全グラフ6頂点, 新規ノードあたり5辺, 有向, シード42
DirectedGraph baGraph = BA.generate("BA", 10000, 6, 5, true, 42L);

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
DirectedGraph g = BA.generate("BA", 1000, 5, 3, true, 42L);

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

### Barabási-Albert モデル（`BA`）

優先的選択（preferential attachment）に基づくスケールフリーネットワークを生成します。

- **有向モード**: 新規ノードは既存ノードの**入次数**に比例する確率で選択
- **無向モード**: 新規ノードは既存ノードの**次数**に比例する確率で選択

```java
// 有向 BA モデル
DirectedGraph gDirected = BA.generate("BA", 10000, 6, 5, true, 42L);

// 無向 BA モデル
DirectedGraph gUndirected = BA.generate("BA_undirected", 10000, 6, 5, false, 42L);
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

---

## 出力ディレクトリ構造

SAR シミュレーションの結果は以下の構造で保存されます:

```
out/fastsar/{path}/{networkType}/threshold={T}/N={N}/{paramKey}={value}/results_{batch_idx}.csv
```

例:
```
out/fastsar/valious_T/DirectedCMOutPow/threshold=4/N=100000/kOutMin=5/results_00.csv
```
