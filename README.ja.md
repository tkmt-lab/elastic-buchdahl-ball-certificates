# 相対論的弾性球の区間証明書

版 **1.0.0 — 2026-11-09**。DOI: —。[論文](paper.pdf)。

Alho・Natário・Pani・Raposo の一つの材料と中心値 δc=8 に対して、中心から有限の表面まで正則な静的弾性球の存在を区間演算で認証する。内部全域で非負圧・優勢エネルギー条件（DEC）・主方向に沿う五つの二乗音速 0≤c²≤1 を満たす。
表面では **0.91089951 < 2M/R < 0.91180552** であり、8/9 を超える。
この材料族で原著者らが示した数値的な超過のうち、一つの異方的弾性球を全域で厳密に確かめる証明書である。Buchdahl の完全流体の仮定とは異なる模型を扱う。

Python 3.10以上、macOS/Linux。展開した一式のルートで実行する：
```sh
python3 -I -B verify.py
```
最後のJSONの期待値： `"status": "PASS"`, `"full_replay": true`。
既定の上限：Aは3300秒、全体は3600秒。[REPRODUCE.md](REPRODUCE.md)を参照。

## 検証と数学的な内容

検査は二つの目録を照合し、証明書を一時ディレクトリに写して、Bの1,592本の管とAの465本の管・4,185成分の微分包囲を全件再計算する。厳密な結果と論文の境界を比較し、すべて通ったときだけ最後に PASS とする。二つの証明書は異なる継続方式で構成され、配布用の再生には共通の整数区間エンジンを使う。Python標準ライブラリでオフライン実行できる。既測定の機械では約8分で、遅い機械ではさらに時間がかかる。上限と失敗時の扱いは [REPRODUCE.md](REPRODUCE.md) に示す。

`python3 -I -B negative_controls.py` は9件の拒否検査。全再現の確認には冒頭のコマンドも実行する。

固定した規格化で半径は 0.01258970 < R < 0.01259407。B は中心の補題の終点包囲を最初の初期箱へ含め、開始時刻も x₀ と一致することを検査する。A も全初期箱からの解を Picard の自己包含で管に閉じ込め、中心から同じ解へ接続する。論文 §3.2 の式(7)の近くに、Bの最初の初期箱の正確な二進端点6個と、実際の中心誤差定数Bを示す。本文の係数再帰と有理数比較だけで包含を照合できる。論文 §3.3 は非自律流の c₁₃(τ,v) と、その時刻までの変分流 V(t) を定義する。

## 範囲

貢献は一つの材料と中心値での全域の厳密認証に限る。基準剪断係数は μ=−60。材料の物理的妥当性、任意方向の音速、安定性、一般の材料についての結論は含まない。現象そのものや厳密認証の先取権は主張しない。Lean は命題の型検査（L1）のみで、定理の形式証明ではない。人間による外部査読は受けていない。下記の Merkushev の計算はこちらで再認証していない。静的AIレビューの範囲と対応履歴は [VERIFICATION-RECORD.md](VERIFICATION-RECORD.md) に示す。

## 対応する版の取得と接続条件

版の取得先は [v1.0.0](https://github.com/tkmt-lab/elastic-buchdahl-ball-certificates/releases/tag/v1.0.0)。

論文が固定する数値入力・検査器の目録 `CERTIFICATE-SHA256SUMS` の SHA-256 は `431b5dc0ad9a67b3f17388b0560fbf6e1556431801ddd0ecc0e51d0e82fd15fa`。macOS の `shasum -a 256 CERTIFICATE-SHA256SUMS` または Linux の `sha256sum CERTIFICATE-SHA256SUMS` で論文記載値と照合する。`verify.py` はこの目録と一式全体の `MANIFEST.sha256` を両方検査する。

## 収録内容

- `certificates/`：二つの証明書、整数区間エンジンと検査器。
- `spec.json`、`FORMAT.md`：固定した材料と証明書の形式。
- `paper.tex`、`paper.pdf`：数学的な議論、定量的な境界、中心と初期箱の接続、Aの管と剰余の定義。
- `MainStatement.lean`、`lean-check.json`：存在命題とL1型検査の記録。
- `CERTIFICATE-SHA256SUMS`、`MANIFEST.sha256`：論文が固定する数値入力・検査器の目録と、全体の目録。
- `REVIEW-RESPONSE.md`、`VERIFICATION-RECORD.md`、`provenance.json`：指摘ごとの対応、検証の根拠と入力の来歴。

## 関連研究と謝意

別の構成則・中心値1/2について、[Merkushev (2026)](https://doi.org/10.5281/zenodo.20836296) は証明書付き静的球を報告する。M/R は約0.11749で、今回のACS材料・高コンパクトネス球とは異なる。

材料・数値例・相対論的弾性論の枠組みを提供した原著者らに感謝する。出典と正確なパラメータ対応は [SOURCES.md](SOURCES.md) に示す。

## ライセンスと連絡先

コードは Apache-2.0、本文・文書・自作データは CC BY 4.0。[LICENSE](LICENSE) を参照。引用情報は [CITATION.cff](CITATION.cff) にある。

滝本竜彦（Tatsuhiko Takimoto）— 研究用連絡先：[research@tatsuhikotakimoto.com](mailto:research@tatsuhikotakimoto.com)。
