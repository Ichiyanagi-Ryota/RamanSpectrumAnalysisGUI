import sys
import os
import numpy as np
import matplotlib.pyplot as plt

from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QListWidget,
    QGroupBox,
    QFileDialog,
    QComboBox,
    QTextEdit,
    QMessageBox,
    QStyle,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.special import expit
from scipy.signal import savgol_filter, find_peaks

# AsLS
def baseline_als(y, lam = 1e7, p = 0.0001, niter = 20):
    L = len(y)

    #2階差分行列
    D = sparse.diags(
        [1.0, -2.0, 1.0],
        offsets = [0, 1, 2],
        shape = (L - 2, L),
        dtype = float
        )

    #重み
    w = np.ones(L, dtype = float)
    for i in range(niter):
        W = sparse.diags(w, offsets = 0, shape = (L, L), dtype = float)
        
        Z = W + lam * D.T @ D

        #計算効率のためCSC形式に変換
        Z = Z.tocsc()

        z = spsolve(Z, w * y)

        #重みの更新
        w = np.where(y > z, p, 1 -p)
    return z

# arPLS
def baseline_arpls(y, lam = 1e7, ratio = 1e-6, niter = 50):
    L = len(y)

    # 2段差分行列
    D = sparse.diags([1.0, -2.0, 1.0], offsets = [0, 1, 2], shape = (L - 2, L), dtype = float)

    # 初期重み
    w = np.ones(L, dtype = float)
    for i in range(niter):
        # 重み行列
        W = sparse.diags(w, offsets = 0, shape = (L, L), dtype = float)
        # ペナルティ付き最小二乗
        Z = W + lam * D.T @ D

        # 計算効率のためのCSC形式
        Z = Z.tocsc()

        # ベースラインの計算
        z = spsolve(Z, w * y)

        # 残差
        d = y - z

        # ベースラインより下側にあるデータ
        dn = d[d < 0]

        # 下側データが存在しない場合
        if len(dn) ==0:
            break

        # 平均と標準偏差
        m = np.mean(dn)
        s = np.std(dn)

        # 標準偏差が０の場合
        if s ==0:
            break

        # 重みの更新
        w_new = expit(-2 * (d - (2 * s - m)) / s)

        # 重みの変化量
        diff = np.linalg.norm(w_new - w) / np.linalg.norm(w)

        # 重みを更新
        w = w_new

        # 収束判定
        if diff < ratio:
            break

    return z

class RamanSpectrumAnalysisGUI(QMainWindow):

    # ウィンドウ設定 =============================================
    def __init__(self):
        super().__init__()

        self.setWindowTitle("RamanSpectrumAnalysisGUI Ver.1.0.0")
        self.resize(1200, 800)

        # メインウィジェット
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.setWindowIcon(QIcon("RamanSpectrumAnalysis.ico"))

        # メインレイアウト
        main_layout = QVBoxLayout()
        central_widget.setLayout(main_layout)

        # 上部タイトル
        title_label = QLabel("Raman Spectrum Analysis GUI")
        title_label.setStyleSheet("font-size: 24px; font-weight: bold;")

        main_layout.addWidget(title_label)

        # 上段パネル 操作 ------------------------------
        upper_layout = QHBoxLayout()

        # 左 操作パネル
        control_panel = QVBoxLayout()

        # DATA グループボックス
        data_group = QGroupBox("DATA")
        data_layout = QVBoxLayout()

        self.file_list = QListWidget()
        self.file_list.itemClicked.connect(self.select_file)

        # Add File と Remove ボタン
        self.add_button = QPushButton(" Add Files")
        self.add_button.setIcon(
            self.style().standardIcon(
                QStyle.StandardPixmap.SP_FileDialogNewFolder
            )
        )
        self.remove_button = QPushButton(" Discard")
        self.remove_button.setIcon(
                self.style().standardIcon(
                QStyle.StandardPixmap.SP_TrashIcon
            )
        )
        self.save_button = QPushButton(" Save Results")
        self.save_button.setIcon(
            self.style().standardIcon(
                QStyle.StandardPixmap.SP_DialogSaveButton
            )
        )
        self.save_all_button = QPushButton(" Save All Results")
        self.save_all_button.setIcon(
            self.style().standardIcon(
                QStyle.StandardPixmap.SP_DialogSaveButton
            )
        )

        # ボタンと機能を接続
        self.add_button.clicked.connect(self.add_files)
        self.remove_button.clicked.connect(self.remove_files)
        self.save_button.clicked.connect(self.save_results)
        self.save_all_button.clicked.connect(self.save_all_results)

        data_layout.addWidget(self.file_list)
        data_layout.addWidget(self.add_button)
        data_layout.addWidget(self.remove_button)
        data_layout.addWidget(self.save_button)
        data_layout.addWidget(self.save_all_button)

        data_group.setLayout(data_layout)

        # ANALYSIS グループボックス
        analysis_group = QGroupBox("ANALYSIS")
        analysis_layout = QVBoxLayout()

        baseline_label = QLabel("Baseline Method")

        # ベースライン方法の選択
        self.baseline_combo = QComboBox()
        self.baseline_combo.addItems([
            "arPLS",
            "AsLS"
        ])

        # 全データへの適用
        self.apply_button = QPushButton("APPLY TO ALL")
        self.apply_button.setIcon(
            self.style().standardIcon(
                QStyle.StandardPixmap.SP_MediaPlay
            )
        )

        self.apply_button.clicked.connect(self.apply_baseline_to_all)

        analysis_layout.addWidget(baseline_label)
        analysis_layout.addWidget(self.baseline_combo)
        analysis_layout.addWidget(self.apply_button)

        analysis_group.setLayout(analysis_layout)

        # 左側パネルに追加
        control_panel.addWidget(data_group)
        control_panel.addWidget(analysis_group)
        control_panel.addStretch()

        # 右 グラフ表示パネル ----------------------------
        graph_group = QGroupBox("SPECTRUM")
        graph_layout = QVBoxLayout()

        # Matplotlib Figure
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)

        # グラフ用
        self.ax = self.figure.add_subplot(111)

        # 解析結果の保存
        self.analysis_results = {}

        # 初期表示
        self.ax.set_xlabel("Raman shift [cm⁻¹]")
        self.ax.set_ylabel("Intensity")
        self.ax.set_title("Raman Spectrum")
        self.ax.grid(True)                  

        graph_layout.addWidget(self.canvas)
        graph_group.setLayout(graph_layout)

        # 左右比率
        upper_layout.addLayout(control_panel, 1)
        upper_layout.addWidget(graph_group, 3)

        main_layout.addLayout(upper_layout)

        # 下段パネル 解析結果 ----------------------------
        result_group = QGroupBox("RESULTS")
        result_layout = QVBoxLayout()

        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)

        result_layout.addWidget(self.result_text)
        result_group.setLayout(result_layout)
        main_layout.addWidget(result_group)

    # ファイル追加 ==============================================
    def add_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Ramanデータを選択してください",
            "",
            "Text files (*.txt);;All files (*)"
        )

        if not file_paths:
            return

        for file_path in file_paths:

            file_name = os.path.basename(file_path)

            # DATA欄に追加
            self.file_list.addItem(file_name)

            # 実際のファイルパスを保存
            list_item = self.file_list.item(self.file_list.count() - 1)
            list_item.setData(Qt.UserRole, file_path)

    # ファイル削除 ==============================================
    def remove_files(self):
        selected_items = self.file_list.selectedItems()

        for item in selected_items:
            row = self.file_list.row(item)
            self.file_list.takeItem(row)

    # ファイルの選択 =============================================
    def select_file(self, item):

        # 実際のファイルパスを取得
        file_path = item.data(Qt.UserRole)

        print("選択されたファイル：")
        print(file_path)

        try:
            # すでに解析済みの場合
            if file_path in self.analysis_results:
                result = self.analysis_results[file_path]

                x = result["x"]
                y = result["raw"]
                baseline = result["baseline"]
                corrected = result["corrected"]
                smoothed = result["smoothed"]

                peak_positions = result["peak_positions"]
                peak_intensities = result["peak_intensities"]

                method = result["method"]

                # グラフをクリア
                self.ax.clear()

                # グラフに描画
                # Raw Spectrum
                self.ax.plot(x, y, label = "Raw Spectrum")
                # Baseline
                self.ax.plot(x, baseline, label = "Baseline")
                # Corrected Spectrum
                self.ax.plot(x, corrected, label = "Corrected Spectrum")
                # Smoothed Spectrum
                self.ax.plot(x, smoothed, label = "Smoothed Spectrum")
                # Detected Peaks
                self.ax.scatter(peak_positions, peak_intensities, marker = "x", s = 60, label = "Detected Peaks")
                
                # 軸ラベル
                self.ax.set_xlabel("Raman Shift [cm⁻¹]")
                self.ax.set_ylabel("Intensity")

                self.ax.set_title(f"{item.text()} - {method}")

                self.ax.legend()
                self.ax.grid(True)

                self.canvas.draw()

                # RESULTSを更新
                self.display_analysis_result(file_path)

            # 未解析の場合
            else:
                # txtファイルを読み込む
                data = np.loadtxt(file_path, dtype = float)

                x = data[:, 0] # Raman Shift
                y = data[:, 1] # Intensity

                # グラフをクリア
                self.ax.clear()

                # スペクトルを描画
                self.ax.plot(x, y)

                # 軸ラベル
                self.ax.set_xlabel("Raman Shift [cm⁻¹]")
                self.ax.set_ylabel("Intensity")

                self.ax.set_title(item.text())

                self.ax.grid(True)

                # グラフを更新
                self.canvas.draw()

                self.result_text.setPlainText(
                    f"Selected File : {item.text()}\n\n"
                    "Analysis has not been performed yet."
                )

        except Exception as e:
            print("ファイル読み込みエラー")
            print(e)

    # 全スペクトルの処理
    def apply_baseline_to_all(self):
        method = self.baseline_combo.currentText()

        if self.file_list.count() == 0:
            print("ファイルがありません")
            return

        print("====================")
        print("Baseline correction started")
        print("Method : ", method)
        print("====================")

        # 解析結果を一旦クリア
        self.analysis_results.clear()

        # 全ファイルを順番に処理
        for i in range(self.file_list.count()):
            item = self.file_list.item(i)

            # ファイルパスを取得
            file_path = item.data(Qt.UserRole)

            print("解析中 : ", item.text())

            try:
                # txtファイル読み込み
                data = np.loadtxt(file_path, dtype = float)

                x = data[:, 0]
                y = data[:, 1]

                # Baseline計算
                if method == "arPLS":
                    baseline = baseline_arpls(y)
                elif method == "AsLS":
                    baseline = baseline_als(y)
                else:
                    print("未検地のベースライン方法です。")
                    continue

                # ベースライン補正
                corrected = y - baseline

                # ピーク検出
                smoothed, peaks, properties = self.detect_peaks(x, corrected)

                peak_positions, peak_intensities = self.get_peak_intensities(
                    x,
                    corrected,
                    peaks
                )

                peak_names = self.assign_peak_names(peak_positions)

                peak_ratios = self.calculate_peak_ratios(
                    peak_names,
                    peak_intensities
                )

                # 結果を保存
                self.analysis_results[file_path] = {
                    "filename": item.text(),
                    "x": x,
                    "raw": y,
                    "baseline": baseline,
                    "corrected": corrected,
                    "method": method,
                    "smoothed": smoothed,
                    "peaks": peaks,
                    "peak_positions": peak_positions,
                    "peak_intensities": peak_intensities,
                    "peak_names": peak_names,
                    "peak_ratios": peak_ratios
                }
                print("完了 : ", item.text())

                print("検出されたピーク")
                for peak_x, peak_intensity, peak_name in zip(
                    peak_positions,
                    peak_intensities,
                    peak_names
                ):
                    print(
                        f"Raman Shift : {peak_x:.1f} cm⁻¹,"
                        f"Intensity : {peak_intensity:.2f}",
                        f"({peak_name})"
                    )

                print("Peak Ratios")
                for ratio_name, ratio_value in peak_ratios.items():
                    if ratio_value is not None:
                        print(f"{ratio_name} : {ratio_value:.4f}")
                    else:
                        print(f"{ratio_name} : N/A")

            except Exception as e:
                print("Baseline correction error:")
                print(item.text())
                print(e)

        print("====================")
        print("Baseline correction finished")
        print("Processed files : ", len(self.analysis_results))
        print("====================")

        # 現在選択されているファイルをグラフ表示
        current_item = self.file_list.currentItem()

        if current_item is not None:
            self.select_file(current_item)

    # ピーク検出
    def detect_peaks(self,x, corrected):
        # Savitzky-Golay smoothing
        smoothed = savgol_filter(
            corrected,
            window_length = 31,
            polyorder = 3
        )

        # Peak detection
        dx = np.median(np.abs(np.diff(x)))

        min_distance_cm = 30

        distance_points = max(
            1,
            int(round(min_distance_cm / dx))
        )

        peaks, properties = find_peaks(
            smoothed,
            prominence=10,
            distance=distance_points
        )
        return smoothed, peaks, properties

    # ピーク強度の特定
    def get_peak_intensities(self, x, corrected, peaks):
        peak_positions = []
        peak_intensities = []

        search_range = 20

        for peak in peaks:

            center_x = x[peak]

            mask = (
                (x >= center_x - search_range) &
                (x <= center_x + search_range)
            )

            # 検索範囲内にデータがない場合
            if not np.any(mask):
                continue

            local_x = x[mask]
            local_y = corrected[mask]

            max_index = np.argmax(local_y)

            peak_x = local_x[max_index]
            peak_intensity = local_y[max_index]

            peak_positions.append(peak_x)
            peak_intensities.append(peak_intensity)

        return(
            np.array(peak_positions),
            np.array(peak_intensities)
        )

    # 主要ピークの命名
    def assign_peak_names(self, peak_positions):
        peak_targets = {
            "Si": 520,
            "D band": 1350,
            "G band": 1580,
            "2D band": 2700
        }

        peak_tolerance = 30

        peak_names = []

        for x_peak in peak_positions:

            assigned_name = "Other"

            for name, target in peak_targets.items():

                if abs(x_peak - target) <= peak_tolerance:

                    assigned_name = name
                    break

            peak_names.append(assigned_name)

        return peak_names

    # ピーク比計算
    def calculate_peak_ratios(self, peak_names, peak_intensities):
        ratio_settings = [
            ("D/G", "D band", "G band"),
            ("2D/G", "2D band", "G band"),
            ("G/Si", "G band", "Si")
            ]

        intensity_dict = {}

        for i in range(len(peak_names)):

            name = peak_names[i]
            intensity = peak_intensities[i]

            if name not in intensity_dict:
                intensity_dict[name] = intensity

        ratio_results = {}

        for ratio_name, numerator_name, denominator_name in ratio_settings:

            if (
                numerator_name in intensity_dict
                and denominator_name in intensity_dict
            ):

                numerator = intensity_dict[numerator_name]
                denominator = intensity_dict[denominator_name]

                if denominator != 0:

                    ratio_results[ratio_name] = (
                        numerator / denominator
                    )

                else:

                    ratio_results[ratio_name] = None

            else:

                ratio_results[ratio_name] = None

        return ratio_results

    # 解析結果のディスプレイ表示
    def display_analysis_result(self, file_path):
        if file_path not in self.analysis_results:
            return
        result = self.analysis_results[file_path]

        filename = result["filename"]
        method = result["method"]

        peak_positions = result["peak_positions"]
        peak_intensities = result["peak_intensities"]
        peak_names = result["peak_names"]
        peak_ratios = result["peak_ratios"]

        text = ""

        text += f"Filename : "
        text += f"{filename}\n"

        text += f"Baseline Method : "
        text += f"{method}\n\n"

        text += "Peak Information =========================\n"
        text += (
            f"{'Raman Shift':<15}"
            f"{'Intensity':<15}"
            f"{'Peak':>8}\n"
        )
        text += "-" * 50 + "\n"

        for peak_name, peak_x, peak_intensity in zip(
            peak_names,
            peak_positions,
            peak_intensities
        ):
            text += (
                f"{peak_x:<20.1f}"
                f"{peak_intensity:<20.2f}"
                f"{peak_name:<12}\n"
            )

        text += "\n"
        text += "Peak Ratios (Graphene) ===================\n"

        for ratio_name, ratio_value in peak_ratios.items():
            if ratio_value is not None:
                text += f"{ratio_name:<10}{ratio_value:.4f}\n"
            else:
                text += f"{ratio_name:<10}N/A\n"

        self.result_text.setPlainText(text)

    # データの保存
    def save_results(self):
        current_item = self.file_list.currentItem()

        if current_item is None:
            QMessageBox.warning(
                self,
                "Save Results",
                "Please select a file."
            )
            return

        file_path = current_item.data(Qt.UserRole)

        if file_path not in self.analysis_results:
            QMessageBox.warning(
                self,
                "Save Results",
                "This file has not been analyzed yet."
            )
            return

        result = self.analysis_results[file_path]

        x = result["x"]
        corrected = result["corrected"]
        smoothed = result["smoothed"]

        peak_positions = result["peak_positions"]
        peak_intensities = result["peak_intensities"]
        peak_names = result["peak_names"]
        peak_ratios = result["peak_ratios"]

        filename = result["filename"]
        method = result["method"]

        # 元フォルダと同じフォルダに保存
        base_name = os.path.splitext(file_path)[0]

        corrected_output = base_name + "_corrected.txt"
        tiff_output = base_name + "_analysis.tiff"

        # Corrected Spectrumのtxt保存
        data_to_save = np.column_stack((x, corrected))

        np.savetxt(
            corrected_output,
            data_to_save,
            fmt = "%.6f",
            delimiter = "\t",
            header = "Raman Shift\tCorrected Intensity",
            comments = ""
        )

        # TIFF保存用グラフの作成
        fig = plt.figure(figsize = (18,9))

        gs = GridSpec(
            1,
            2,
            width_ratios = [5, 1],
            wspace = 0.01,
            figure = fig
            )

        ax = fig.add_subplot(gs[0])
        info_ax = fig.add_subplot(gs[1])

        info_ax.axis("off")

        # スペクトル
        ax.plot(
            x,
            corrected,
            label="Corrected Spectrum"
        )

        ax.plot(
            x,
            smoothed,
            label="Smoothed Spectrum"
        )

        # ピーク位置
        ax.scatter(
            peak_positions,
            peak_intensities,
            marker="x",
            label="Detected Peaks"
        )

        ax.set_xlabel("Raman Shift [cm$^{-1}$]")

        ax.set_ylabel("Intensity [a.u.]")

        ax.set_title(filename)

        ax.legend()
        ax.grid(True)

        # 右側の情報欄
        info_ax.axis("off")
        info_text = ""

        info_text += ("Original File : " + f"{filename}\n\n")

        info_text += ("Baseline Method : " + f"{method}\n\n")

        info_text += "Peak Information --------------------\n\n"

        for i in range(len(peak_positions)):

            info_text += (
                f"{peak_positions[i] : .2f} cm⁻¹  "
                f"{peak_intensities[i] : .2f}  "
                f"({peak_names[i]})\n"
            )

        info_text += "\nPeak Ratios (Graphene) --------------------\n\n"

        for ratio_name, ratio_value in peak_ratios.items():

            if ratio_value is not None:
                info_text += (f"{ratio_name} : " f"{ratio_value:.4f}\n")
            else:
                info_text += (f"{ratio_name} : N/A\n")

        info_ax.text(
            0,
            1,
            info_text,
            transform=info_ax.transAxes,
            verticalalignment="top",
            fontsize=10
        )

        # TIFF保存
        plt.savefig(
            tiff_output,
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(fig)

        QMessageBox.information(
            self,
            "Save Results",
            "Results have been saved successfully.\n\n"
            f"Corrected Spectrum:\n{corrected_output}\n\n"
            f"TIFF:\n{tiff_output}"
        )

    # 全データを一括保存
    def save_all_results(self):
        if not self.analysis_results:
            QMessageBox.earning(
                self,
                "Save All Results",
                "No analyzed results are available."
            )
            return

        saved_count = 0

        for file_path, result in self.analysis_results.items():
            x = result["x"]
            corrected = result["corrected"]
            smoothed = result["smoothed"]

            peak_positions = result["peak_positions"]
            peak_intensities = result["peak_intensities"]
            peak_names = result["peak_names"]
            peak_ratios = result["peak_ratios"]

            filename = result["filename"]
            method = result["method"]

            base_name = os.path.splitext(file_path)[0]

            corrected_output = base_name + "_corrected.txt"
            tiff_output = base_name + "_analysis.tiff"

            data_to_save = np.column_stack((x, corrected))
    
            np.savetxt(
                corrected_output,
                data_to_save,
                fmt = "%.6f",
                delimiter = "\t",
                header = "Raman Shift\tCorrected Intensity",
                comments = ""
            )
    
            # TIFF保存用グラフの作成
            fig = plt.figure(figsize = (18,9))
    
            gs = GridSpec(
                1,
                2,
                width_ratios = [5, 1],
                wspace = 0.01,
                figure = fig
                )
    
            ax = fig.add_subplot(gs[0])
            info_ax = fig.add_subplot(gs[1])
    
            info_ax.axis("off")
    
            # スペクトル
            ax.plot(
                x,
                corrected,
                label="Corrected Spectrum"
            )
    
            ax.plot(
                x,
                smoothed,
                label="Smoothed Spectrum"
            )
    
            # ピーク位置
            ax.scatter(
                peak_positions,
                peak_intensities,
                marker="x",
                label="Detected Peaks"
            )
    
            ax.set_xlabel("Raman Shift [cm$^{-1}$]")
    
            ax.set_ylabel("Intensity [a.u.]")
    
            ax.set_title(filename)
    
            ax.legend()
            ax.grid(True)
    
            # 右側の情報欄
            info_ax.axis("off")
            info_text = ""
    
            info_text += ("Original File : " + f"{filename}\n\n")
    
            info_text += ("Baseline Method : " + f"{method}\n\n")
    
            info_text += "Peak Information --------------------\n\n"
    
            for i in range(len(peak_positions)):
    
                info_text += (
                    f"{peak_positions[i] : .2f} cm⁻¹  "
                    f"{peak_intensities[i] : .2f}  "
                    f"({peak_names[i]})\n"
                )
    
            info_text += "\nPeak Ratios (Graphene) --------------------\n\n"
    
            for ratio_name, ratio_value in peak_ratios.items():
    
                if ratio_value is not None:
                    info_text += (f"{ratio_name} : " f"{ratio_value:.4f}\n")
                else:
                    info_text += (f"{ratio_name} : N/A\n")
    
            info_ax.text(
                0,
                1,
                info_text,
                transform=info_ax.transAxes,
                verticalalignment="top",
                fontsize=10
            )
    
            # TIFF保存
            plt.savefig(
                tiff_output,
                dpi=300,
                bbox_inches="tight"
            )
    
            plt.close(fig)

            saved_count += 1
    
        QMessageBox.information(
            self,
            "Save All Results",
            f"{saved_count} file(s) have been saved succesfully"
        )

def main():
    app = QApplication(sys.argv)

    app.setWindowIcon(QIcon("RamanSpectrumAnalysis.ico"))

    window = RamanSpectrumAnalysisGUI()
    window.show()

    sys.exit(app.exec())

if __name__ == "__main__":
    main()