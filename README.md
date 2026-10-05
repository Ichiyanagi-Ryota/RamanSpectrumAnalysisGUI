RamanSpectrumAnalysisGUI
Version 1.0.0

【概要】
RamanSpectrumAnalysisGUIは、Ramanスペクトルの
読み込み、Baseline補正、Peak Detection、
Peak Intensity取得および解析結果の保存を行うGUIアプリケーションです。

【起動方法】
RamanSpectrumAnalysisGUI.exeをダブルクリックしてください。

【基本操作】

1. Add Files
   Ramanスペクトルのtxtファイルを選択します。

2. DATA欄から解析対象のファイルを選択します。

3. Baseline Methodを選択します。
   - arPLS
   - AsLS

4. APPLY TO ALLを押します。

5. 解析結果を確認します。

【結果の保存】

Save Results
現在選択しているファイルの結果を保存します。

Save All Results
解析済みの全ファイルの結果を保存します。

保存されるファイル：
・*_corrected.txt
・*_analysis.tiff

【注意】
入力データは2列の数値データを想定しています。

1列目：Raman Shift
2列目：Intensity