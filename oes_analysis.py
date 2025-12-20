#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OES (Optical Emission Spectroscopy) 데이터 분석 프로그램
수소 플라즈마 발머 계열 파장 분석 및 시각화
"""

import sys
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QPushButton, QLabel, QDoubleSpinBox,
                             QFileDialog, QMessageBox, QSplitter)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


class OESAnalyzer(QMainWindow):
    """OES 데이터 분석 메인 윈도우"""

    # 발머 계열 파장 정의 (nm)
    BALMER_WAVELENGTHS = {
        'Hα': 656.28,
        'Hβ': 486.13,
        'Hγ': 434.05
    }

    def __init__(self):
        super().__init__()
        self.data = None  # 로딩된 데이터프레임
        self.run_times = None  # 시간 배열
        self.wavelengths = None  # 파장 배열
        self.current_time = 0.0  # 현재 선택된 시간
        self.window_size = 2.0  # 가우시안 윈도우 크기 (기본값 ±2nm)
        self.time_line = None  # 창2의 시간 표시 선

        self.init_ui()

    def init_ui(self):
        """GUI 초기화"""
        self.setWindowTitle('OES 데이터 분석 프로그램')
        self.setGeometry(100, 100, 1400, 900)

        # 중앙 위젯 생성
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # 메인 레이아웃 (수평 분할: 컨트롤 패널 1 : 그래프 영역 4)
        main_layout = QHBoxLayout(central_widget)

        # 왼쪽 컨트롤 패널 생성
        control_panel = self.create_control_panel()

        # 오른쪽 그래프 영역 생성 (상하 분할)
        graph_widget = self.create_graph_area()

        # 수평 분할기로 비율 설정 (1:4)
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(control_panel)
        splitter.addWidget(graph_widget)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 4)

        main_layout.addWidget(splitter)

    def create_control_panel(self):
        """왼쪽 컨트롤 패널 생성"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setAlignment(Qt.AlignTop)

        # 파일 로딩 버튼
        self.load_button = QPushButton('파일 로딩')
        self.load_button.clicked.connect(self.load_file)
        layout.addWidget(self.load_button)

        layout.addSpacing(20)

        # Window 입력
        window_label = QLabel('Window (nm)')
        layout.addWidget(window_label)

        self.window_spinbox = QDoubleSpinBox()
        self.window_spinbox.setRange(0.5, 10.0)
        self.window_spinbox.setSingleStep(0.5)
        self.window_spinbox.setValue(2.0)
        self.window_spinbox.setDecimals(1)
        self.window_spinbox.valueChanged.connect(self.on_window_changed)
        layout.addWidget(self.window_spinbox)

        layout.addSpacing(20)

        # 시간 입력
        time_label = QLabel('Time (sec)')
        layout.addWidget(time_label)

        self.time_spinbox = QDoubleSpinBox()
        self.time_spinbox.setEnabled(False)  # 초기에는 비활성화
        self.time_spinbox.valueChanged.connect(self.on_time_changed)
        layout.addWidget(self.time_spinbox)

        layout.addStretch()

        return panel

    def create_graph_area(self):
        """오른쪽 그래프 영역 생성 (상하 분할)"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)

        # 창1: 스펙트럼 그래프
        self.spectrum_figure = Figure(figsize=(8, 4), facecolor='white')
        self.spectrum_canvas = FigureCanvas(self.spectrum_figure)
        self.spectrum_ax = self.spectrum_figure.add_subplot(111)
        self.spectrum_toolbar = NavigationToolbar(self.spectrum_canvas, widget)

        spectrum_widget = QWidget()
        spectrum_layout = QVBoxLayout(spectrum_widget)
        spectrum_layout.setContentsMargins(0, 0, 0, 0)
        spectrum_layout.addWidget(self.spectrum_toolbar)
        spectrum_layout.addWidget(self.spectrum_canvas)

        # 창2: 시계열 그래프
        self.timeseries_figure = Figure(figsize=(8, 4), facecolor='white')
        self.timeseries_canvas = FigureCanvas(self.timeseries_figure)
        self.timeseries_ax = self.timeseries_figure.add_subplot(111)
        self.timeseries_toolbar = NavigationToolbar(self.timeseries_canvas, widget)

        # 창2에 클릭 이벤트 연결
        self.timeseries_canvas.mpl_connect('button_press_event', self.on_timeseries_click)

        timeseries_widget = QWidget()
        timeseries_layout = QVBoxLayout(timeseries_widget)
        timeseries_layout.setContentsMargins(0, 0, 0, 0)
        timeseries_layout.addWidget(self.timeseries_toolbar)
        timeseries_layout.addWidget(self.timeseries_canvas)

        # 수직 분할기로 1:1 비율 설정
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(spectrum_widget)
        splitter.addWidget(timeseries_widget)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)

        layout.addWidget(splitter)

        # 초기 빈 그래프 설정
        self.setup_empty_graphs()

        return widget

    def setup_empty_graphs(self):
        """초기 빈 그래프 설정"""
        # 창1: 스펙트럼 그래프
        self.spectrum_ax.clear()
        self.spectrum_ax.set_xlabel('Wavelength (nm)', fontsize=10)
        self.spectrum_ax.set_ylabel('Emission Intensity (a.u.)', fontsize=10)
        self.spectrum_ax.set_title('Spectrum', fontsize=10)
        self.spectrum_ax.grid(True, alpha=0.7, linestyle='--', color='lightgray')
        self.spectrum_canvas.draw()

        # 창2: 시계열 그래프
        self.timeseries_ax.clear()
        self.timeseries_ax.set_xlabel('Run Time (sec)', fontsize=10)
        self.timeseries_ax.set_ylabel('Intensity (a.u.)', fontsize=10)
        self.timeseries_ax.set_title('Balmer Series Time Trace', fontsize=10)
        self.timeseries_ax.grid(True, alpha=0.7, linestyle='--', color='lightgray')
        self.timeseries_canvas.draw()

    def load_file(self):
        """파일 로딩"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, '파일 선택', '', 'Data Files (*.dat);;All Files (*)'
        )

        if not file_path:
            return

        try:
            # 데이터 파일 읽기 (탭 구분, 첫 행은 헤더)
            self.data = pd.read_csv(file_path, sep='\t', encoding='utf-8')

            # 파장 배열 생성 (200.0nm ~ 800.0nm, 0.5nm 간격)
            self.wavelengths = np.arange(200.0, 800.5, 0.5)

            # Run Time 추출 (두 번째 열)
            self.run_times = self.data.iloc[:, 1].values

            # 시간 SpinBox 설정
            time_step = np.min(np.diff(self.run_times)) if len(self.run_times) > 1 else 0.1
            self.time_spinbox.setRange(self.run_times.min(), self.run_times.max())
            self.time_spinbox.setSingleStep(time_step)
            self.time_spinbox.setDecimals(2)
            self.time_spinbox.setValue(self.run_times[0])
            self.time_spinbox.setEnabled(True)

            # 첫 번째 시간으로 초기화
            self.current_time = self.run_times[0]

            # 그래프 업데이트
            self.update_spectrum()
            self.update_timeseries()

            QMessageBox.information(self, '성공', '파일이 성공적으로 로딩되었습니다.')

        except Exception as e:
            QMessageBox.critical(self, '오류', f'파일 로딩 실패:\n{str(e)}')

    def gaussian_weighted_average(self, wavelength_center, spectrum_data):
        """
        가우시안 가중 평균 계산

        Parameters:
        -----------
        wavelength_center : float
            중심 파장 (nm)
        spectrum_data : numpy.array
            전체 스펙트럼 데이터 (파장별 강도)

        Returns:
        --------
        float : 가중 평균된 강도 값
        """
        # 시그마 계산: σ = Window / 2
        sigma = self.window_size / 2.0

        # 윈도우 범위 내의 파장 인덱스 찾기
        wavelength_min = wavelength_center - self.window_size
        wavelength_max = wavelength_center + self.window_size

        # 파장 범위 내의 마스크 생성
        mask = (self.wavelengths >= wavelength_min) & (self.wavelengths <= wavelength_max)
        wavelengths_in_window = self.wavelengths[mask]
        intensities_in_window = spectrum_data[mask]

        if len(wavelengths_in_window) == 0:
            return 0.0

        # 가우시안 가중치 계산: w(λ) = exp(-((λ - λ₀)² / (2σ²)))
        weights = np.exp(-((wavelengths_in_window - wavelength_center) ** 2) / (2 * sigma ** 2))

        # 가중 평균: Σ(w(λ) × I(λ)) / Σ(w(λ))
        weighted_avg = np.sum(weights * intensities_in_window) / np.sum(weights)

        return weighted_avg

    def get_spectrum_at_time(self, time):
        """
        특정 시간의 스펙트럼 데이터 추출

        Parameters:
        -----------
        time : float
            추출할 시간 (초)

        Returns:
        --------
        numpy.array : 파장별 강도 배열
        """
        # 가장 가까운 시간의 인덱스 찾기
        idx = np.argmin(np.abs(self.run_times - time))

        # 해당 시간의 스펙트럼 데이터 추출 (3번째 열부터 끝까지)
        spectrum = self.data.iloc[idx, 2:].values.astype(float)

        return spectrum

    def update_spectrum(self):
        """창1: 스펙트럼 그래프 업데이트"""
        if self.data is None:
            return

        # 현재 시간의 스펙트럼 데이터 가져오기
        spectrum = self.get_spectrum_at_time(self.current_time)

        # 그래프 업데이트
        self.spectrum_ax.clear()
        self.spectrum_ax.plot(self.wavelengths, spectrum, color='#1f77b4', linewidth=1.5)
        self.spectrum_ax.set_xlabel('Wavelength (nm)', fontsize=10)
        self.spectrum_ax.set_ylabel('Emission Intensity (a.u.)', fontsize=10)
        self.spectrum_ax.set_title(f'Spectrum at t = {self.current_time:.1f} s', fontsize=10)
        self.spectrum_ax.grid(True, alpha=0.7, linestyle='--', color='lightgray')
        self.spectrum_figure.tight_layout()
        self.spectrum_canvas.draw()

    def update_timeseries(self):
        """창2: 시계열 그래프 업데이트"""
        if self.data is None:
            return

        # 각 발머 계열 파장에 대한 시계열 데이터 계산
        balmer_timeseries = {name: [] for name in self.BALMER_WAVELENGTHS.keys()}

        for i in range(len(self.run_times)):
            spectrum = self.data.iloc[i, 2:].values.astype(float)

            for name, wavelength in self.BALMER_WAVELENGTHS.items():
                intensity = self.gaussian_weighted_average(wavelength, spectrum)
                balmer_timeseries[name].append(intensity)

        # 그래프 업데이트
        self.timeseries_ax.clear()

        # 3개 라인 그리기
        colors = ['C0', 'C1', 'C2']  # Matplotlib 기본 색상 순환
        for i, (name, wavelength) in enumerate(self.BALMER_WAVELENGTHS.items()):
            self.timeseries_ax.plot(
                self.run_times,
                balmer_timeseries[name],
                label=f'{name} ({wavelength} nm)',
                color=colors[i],
                linewidth=1.5
            )

        # 현재 선택된 시간 표시 (빨간 수직 점선)
        self.time_line = self.timeseries_ax.axvline(
            x=self.current_time,
            color='red',
            linestyle='--',
            linewidth=1.5
        )

        self.timeseries_ax.set_xlabel('Run Time (sec)', fontsize=10)
        self.timeseries_ax.set_ylabel('Intensity (a.u.)', fontsize=10)
        self.timeseries_ax.set_title('Balmer Series Time Trace', fontsize=10)
        self.timeseries_ax.legend(loc='best', fontsize=9)
        self.timeseries_ax.grid(True, alpha=0.7, linestyle='--', color='lightgray')
        self.timeseries_figure.tight_layout()
        self.timeseries_canvas.draw()

    def on_time_changed(self, value):
        """시간 SpinBox 값 변경 이벤트"""
        self.current_time = value
        self.update_spectrum()
        self.update_time_line()

    def on_window_changed(self, value):
        """Window SpinBox 값 변경 이벤트"""
        self.window_size = value
        if self.data is not None:
            self.update_timeseries()

    def on_timeseries_click(self, event):
        """창2 그래프 클릭 이벤트"""
        if event.inaxes != self.timeseries_ax:
            return

        if self.data is None:
            return

        # 클릭한 x좌표(시간) 가져오기
        clicked_time = event.xdata

        # 가장 가까운 실제 시간으로 조정
        idx = np.argmin(np.abs(self.run_times - clicked_time))
        selected_time = self.run_times[idx]

        # SpinBox 업데이트 (이것이 자동으로 on_time_changed 호출)
        self.time_spinbox.setValue(selected_time)

    def update_time_line(self):
        """창2의 시간 표시선 업데이트"""
        if self.time_line is not None:
            self.time_line.set_xdata([self.current_time, self.current_time])
            self.timeseries_canvas.draw()


def main():
    """메인 함수"""
    app = QApplication(sys.argv)
    window = OESAnalyzer()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
