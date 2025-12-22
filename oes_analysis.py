#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OES (Optical Emission Spectroscopy) 데이터 분석 프로그램
수소 플라즈마 발머 계열 파장 분석 및 시각화
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import linregress
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QPushButton, QLabel, QDoubleSpinBox,
                             QFileDialog, QMessageBox, QSplitter, QLineEdit, QCheckBox)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class OESAnalyzer(QMainWindow):
    """OES 데이터 분석 메인 윈도우"""

    # 발머 계열 파장 정의 (nm)
    BALMER_WAVELENGTHS = {
        'Hα': 656.28,
        'Hβ': 486.13,
        'Hγ': 434.05
    }

    # NIST 원자 스펙트럼 상수 (Boltzmann Plot Method용)
    BALMER_CONSTANTS = {
        'Hα': {
            'wavelength': 656.28,  # nm
            'g': 18,               # 통계적 가중치 (n=3, 2n²)
            'A': 4.41e7,           # 전이 확률 (s⁻¹)
            'E': 12.09             # 상위 에너지 레벨 (eV)
        },
        'Hβ': {
            'wavelength': 486.13,  # nm
            'g': 32,               # 통계적 가중치 (n=4, 2n²)
            'A': 8.42e6,           # 전이 확률 (s⁻¹)
            'E': 12.75             # 상위 에너지 레벨 (eV)
        },
        'Hγ': {
            'wavelength': 434.05,  # nm
            'g': 50,               # 통계적 가중치 (n=5, 2n²)
            'A': 2.53e6,           # 전이 확률 (s⁻¹)
            'E': 13.05             # 상위 에너지 레벨 (eV)
        }
    }

    # 볼츠만 상수 (eV/K)
    K_B = 8.617333262e-5

    def __init__(self):
        super().__init__()
        self.data = None  # 로딩된 데이터프레임
        self.run_times = None  # 시간 배열
        self.wavelengths = None  # 파장 배열
        self.current_time = 0.0  # 현재 선택된 시간
        self.window_size = 2.0  # 가우시안 윈도우 크기 (기본값 ±2nm)
        self.time_line = None  # 창2의 시간 표시 선
        self.intensity_markers = []  # 창2 클릭 시 표시되는 마커들
        self.intensity_annotations = []  # 창2 클릭 시 표시되는 텍스트들
        self.texc_annotation = None  # 창2에 표시되는 Texc 텍스트
        self.timeseries_ax2 = None  # 창2의 보조 Y축 (Texc용)
        self.balmer_timeseries = {}  # 전체 시계열 Intensity 데이터 캐시
        self.texc_timeseries = []  # 전체 시계열 Texc 데이터
        self.boltzmann_window = None  # 창3: Boltzmann Plot 팝업창

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

        layout.addSpacing(20)

        # 여기 온도 표시
        texc_label = QLabel('여기 온도')
        layout.addWidget(texc_label)

        self.texc_lineedit = QLineEdit()
        self.texc_lineedit.setReadOnly(True)
        self.texc_lineedit.setText('---')
        layout.addWidget(self.texc_lineedit)

        self.r2_label = QLabel('(R² = ---)')
        layout.addWidget(self.r2_label)

        layout.addSpacing(20)

        # Boltzmann Plot 체크박스
        self.boltzmann_checkbox = QCheckBox('Boltzmann Plot')
        self.boltzmann_checkbox.setChecked(False)
        self.boltzmann_checkbox.stateChanged.connect(self.on_boltzmann_checkbox_changed)
        layout.addWidget(self.boltzmann_checkbox)

        layout.addStretch()

        return panel

    def create_graph_area(self):
        """오른쪽 그래프 영역 생성 (상하 분할)"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)

        # 창1: 스펙트럼 그래프 (툴바 제거)
        self.spectrum_figure = Figure(figsize=(8, 4), facecolor='white')
        self.spectrum_canvas = FigureCanvas(self.spectrum_figure)
        self.spectrum_ax = self.spectrum_figure.add_subplot(111)

        spectrum_widget = QWidget()
        spectrum_layout = QVBoxLayout(spectrum_widget)
        spectrum_layout.setContentsMargins(0, 0, 0, 0)
        spectrum_layout.addWidget(self.spectrum_canvas)

        # 창2: 시계열 그래프 (툴바 제거)
        self.timeseries_figure = Figure(figsize=(8, 4), facecolor='white')
        self.timeseries_canvas = FigureCanvas(self.timeseries_figure)
        self.timeseries_ax = self.timeseries_figure.add_subplot(111)

        # 창2에 클릭 이벤트 연결
        self.timeseries_canvas.mpl_connect('button_press_event', self.on_timeseries_click)

        timeseries_widget = QWidget()
        timeseries_layout = QVBoxLayout(timeseries_widget)
        timeseries_layout.setContentsMargins(0, 0, 0, 0)
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

            # 데이터 유효성 검증
            if len(self.run_times) == 0:
                raise ValueError("데이터에 유효한 시간 정보가 없습니다.")

            # 시간 SpinBox 설정 (실제 데이터 범위로 제한)
            min_time = self.run_times.min()
            max_time = self.run_times.max()
            time_step = np.min(np.diff(self.run_times)) if len(self.run_times) > 1 else 0.1

            self.time_spinbox.setMinimum(min_time)
            self.time_spinbox.setMaximum(max_time)
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

    def get_time_index(self, time_value):
        """
        시간 값을 데이터 인덱스로 변환 (범위 검증 포함)

        Parameters:
        -----------
        time_value : float
            변환할 시간 값 (초)

        Returns:
        --------
        int : 유효한 인덱스 (0 ~ len(run_times)-1)
        """
        if self.run_times is None or len(self.run_times) == 0:
            return 0

        # 가장 가까운 시간의 인덱스 찾기
        idx = np.argmin(np.abs(self.run_times - time_value))

        # 범위 검증 (안전장치)
        max_idx = len(self.run_times) - 1
        idx = max(0, min(idx, max_idx))

        return idx

    def calculate_texc(self, intensities):
        """
        Boltzmann Plot Method로 여기 전자 온도 계산

        Parameters:
        -----------
        intensities : dict
            {'Hα': I_alpha, 'Hβ': I_beta, 'Hγ': I_gamma} 형태의 Intensity 딕셔너리

        Returns:
        --------
        tuple : (Texc_eV, R_squared, error_message)
            - Texc_eV: 여기 전자 온도 (eV), 계산 실패 시 None
            - R_squared: 결정 계수, 계산 실패 시 None
            - error_message: 오류 메시지, 성공 시 None
        """
        try:
            # 데이터 준비
            x_data = []  # E_n (eV)
            y_data = []  # ln(I × λ / (g × A))

            for name in ['Hα', 'Hβ', 'Hγ']:
                I = intensities[name]
                const = self.BALMER_CONSTANTS[name]

                # Intensity 유효성 검사
                if I <= 0:
                    return None, None, "계산 불가: 유효하지 않은 Intensity"

                # y = ln(I × λ / (g × A)) 계산
                lambda_nm = const['wavelength']
                g = const['g']
                A = const['A']
                E = const['E']

                y = np.log(I * lambda_nm / (g * A))

                x_data.append(E)
                y_data.append(y)

            # 선형 회귀 (y = slope * x + intercept)
            x_data = np.array(x_data)
            y_data = np.array(y_data)

            # numpy polyfit 사용 (1차 다항식)
            coeffs = np.polyfit(x_data, y_data, 1)
            slope = coeffs[0]
            intercept = coeffs[1]

            # 기울기 검사 (양수면 비물리적)
            if slope >= 0:
                return None, None, "계산 불가: 비물리적 기울기 (양수)"

            # R² 계산
            y_fit = slope * x_data + intercept
            ss_res = np.sum((y_data - y_fit) ** 2)
            ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)

            # Texc 계산: slope = -1 / (k_B × T_exc)
            # T_exc (K) = -1 / (k_B × slope)
            T_exc_K = -1.0 / (self.K_B * slope)

            # eV 단위로 변환: T_exc_eV = k_B × T_exc (K)
            T_exc_eV = self.K_B * T_exc_K

            return T_exc_eV, r_squared, None

        except Exception as e:
            return None, None, f"계산 불가: {str(e)}"

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
        self.balmer_timeseries = {name: [] for name in self.BALMER_WAVELENGTHS.keys()}
        self.texc_timeseries = []

        for i in range(len(self.run_times)):
            spectrum = self.data.iloc[i, 2:].values.astype(float)

            # Intensity 계산
            intensities = {}
            for name, wavelength in self.BALMER_WAVELENGTHS.items():
                intensity = self.gaussian_weighted_average(wavelength, spectrum)
                self.balmer_timeseries[name].append(intensity)
                intensities[name] = intensity

            # Texc 계산
            texc_eV, r2, error_msg = self.calculate_texc(intensities)
            if error_msg is None:
                self.texc_timeseries.append(texc_eV)
            else:
                self.texc_timeseries.append(np.nan)  # 계산 불가 시 NaN

        # 그래프 업데이트
        self.timeseries_ax.clear()

        # 보조 Y축이 있으면 제거
        if self.timeseries_ax2 is not None:
            self.timeseries_ax2.remove()
            self.timeseries_ax2 = None

        # 3개 Intensity 라인 그리기 (좌측 Y축)
        colors = ['C0', 'C1', 'C2']  # Matplotlib 기본 색상 순환
        lines1 = []
        labels1 = []
        for i, (name, wavelength) in enumerate(self.BALMER_WAVELENGTHS.items()):
            line, = self.timeseries_ax.plot(
                self.run_times,
                self.balmer_timeseries[name],
                label=f'{name} ({wavelength} nm)',
                color=colors[i],
                linewidth=1.5
            )
            lines1.append(line)
            labels1.append(f'{name} ({wavelength} nm)')

        # 보조 Y축 생성 (우측)
        self.timeseries_ax2 = self.timeseries_ax.twinx()

        # Texc 라인 그리기 (우측 Y축)
        line2, = self.timeseries_ax2.plot(
            self.run_times,
            self.texc_timeseries,
            label='Texc',
            color='#000000',
            linestyle='-',
            marker='s',
            markersize=4,
            linewidth=1.5
        )

        # 현재 선택된 시간 표시 (빨간 수직 점선)
        self.time_line = self.timeseries_ax.axvline(
            x=self.current_time,
            color='red',
            linestyle='--',
            linewidth=1.5
        )

        # 축 라벨 설정
        self.timeseries_ax.set_xlabel('Run Time (sec)', fontsize=10)
        self.timeseries_ax.set_ylabel('Intensity (a.u.)', fontsize=10)
        self.timeseries_ax2.set_ylabel('Excitation Temperature (eV)', fontsize=10)
        self.timeseries_ax2.tick_params(axis='y', labelcolor='#000000')

        # 범례 통합
        lines = lines1 + [line2]
        labels = labels1 + ['Texc']
        self.timeseries_ax.legend(lines, labels, loc='best', fontsize=9)

        self.timeseries_ax.set_title('Balmer Series Time Trace', fontsize=10)
        self.timeseries_ax.grid(True, alpha=0.7, linestyle='--', color='lightgray')
        self.timeseries_figure.tight_layout()
        self.timeseries_canvas.draw()

        # 초기 Texc 표시 업데이트
        self.update_texc_display()

    def update_texc_display(self):
        """컨트롤 패널의 Texc 표시 업데이트"""
        if self.data is None:
            return

        # 현재 시간의 인덱스 찾기 (범위 검증 포함)
        idx = self.get_time_index(self.current_time)

        # 현재 시간의 Intensity 가져오기
        intensities = {}

        # balmer_timeseries가 비어있으면 직접 계산
        if not self.balmer_timeseries:
            spectrum = self.get_spectrum_at_time(self.current_time)
            for name, wavelength in self.BALMER_WAVELENGTHS.items():
                intensity = self.gaussian_weighted_average(wavelength, spectrum)
                intensities[name] = intensity
        else:
            # 캐시된 데이터 사용 (추가 범위 검증)
            for name in self.BALMER_WAVELENGTHS.keys():
                # balmer_timeseries의 길이 확인
                if idx >= len(self.balmer_timeseries[name]):
                    # 범위 초과 시 마지막 인덱스로 조정
                    idx = len(self.balmer_timeseries[name]) - 1
                intensities[name] = self.balmer_timeseries[name][idx]

        # Texc 계산
        texc_eV, r2, error_msg = self.calculate_texc(intensities)

        if error_msg is None:
            # 성공적으로 계산됨
            self.texc_lineedit.setText(f'{texc_eV:.2f} eV')
            self.r2_label.setText(f'(R² = {r2:.2f})')
        else:
            # 계산 실패
            self.texc_lineedit.setText(error_msg)
            self.r2_label.setText('(R² = ---)')

    def on_time_changed(self, value):
        """시간 SpinBox 값 변경 이벤트"""
        self.current_time = value
        self.update_spectrum()
        self.clear_intensity_markers()  # 마커 초기화
        self.update_time_line()
        self.update_texc_display()
        self.timeseries_canvas.draw()  # 캔버스 업데이트

        # 창3이 열려있으면 업데이트
        self.update_boltzmann_plot()

    def on_window_changed(self, value):
        """Window SpinBox 값 변경 이벤트"""
        self.window_size = value
        if self.data is not None:
            self.update_timeseries()
            # 창3이 열려있으면 업데이트
            self.update_boltzmann_plot()

    def clear_intensity_markers(self):
        """창2의 Intensity 마커와 annotation 제거"""
        # 마커 제거
        for marker in self.intensity_markers:
            marker.remove()
        self.intensity_markers = []

        # Annotation 제거
        for ann in self.intensity_annotations:
            ann.remove()
        self.intensity_annotations = []

        # Texc annotation 제거
        if self.texc_annotation is not None:
            self.texc_annotation.remove()
            self.texc_annotation = None

    def display_intensity_markers(self):
        """창2에 현재 시간의 Intensity 마커와 Texc 표시"""
        if self.data is None or not self.balmer_timeseries:
            return

        # 현재 시간의 인덱스 찾기 (범위 검증 포함)
        idx = self.get_time_index(self.current_time)

        # 각 라인의 Intensity 값 가져오기
        colors = ['C0', 'C1', 'C2']
        intensities = {}

        for i, name in enumerate(['Hα', 'Hβ', 'Hγ']):
            # 인덱스 범위 검증
            if idx >= len(self.balmer_timeseries[name]):
                idx = len(self.balmer_timeseries[name]) - 1
            intensity = self.balmer_timeseries[name][idx]
            intensities[name] = intensity

            # 마커 표시
            marker, = self.timeseries_ax.plot(
                self.current_time,
                intensity,
                marker='o',
                markersize=8,
                color=colors[i],
                markeredgecolor='white',
                markeredgewidth=1.5,
                zorder=10
            )
            self.intensity_markers.append(marker)

            # Annotation 표시
            ann = self.timeseries_ax.annotate(
                f'{intensity:.1f}',
                xy=(self.current_time, intensity),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=9,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
            )
            self.intensity_annotations.append(ann)

        # Texc 계산 및 그래프 내부 표시
        texc_eV, r2, error_msg = self.calculate_texc(intensities)

        if error_msg is None:
            texc_text = f'Texc = {texc_eV:.2f} eV\n(R² = {r2:.2f})'
        else:
            texc_text = error_msg

        # 그래프 우측 상단에 Texc 표시
        self.texc_annotation = self.timeseries_ax.text(
            0.98, 0.98,
            texc_text,
            transform=self.timeseries_ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
        )

        self.timeseries_canvas.draw()

    def on_timeseries_click(self, event):
        """창2 그래프 클릭 이벤트"""
        # 보조 Y축도 포함하여 클릭 가능하도록 설정
        valid_axes = [self.timeseries_ax]
        if self.timeseries_ax2 is not None:
            valid_axes.append(self.timeseries_ax2)

        if event.inaxes not in valid_axes:
            return

        if self.data is None:
            return

        # 클릭한 x좌표(시간) 가져오기
        clicked_time = event.xdata

        if clicked_time is None:
            return

        # 클릭한 시간이 데이터 범위 내인지 검증
        if len(self.run_times) > 0:
            min_time = self.run_times.min()
            max_time = self.run_times.max()
            clicked_time = max(min_time, min(clicked_time, max_time))

        # 가장 가까운 실제 시간으로 조정 (범위 검증 포함)
        idx = self.get_time_index(clicked_time)
        selected_time = self.run_times[idx]

        # 이전 마커와 annotation 제거
        self.clear_intensity_markers()

        # 새로운 시간 설정
        self.current_time = selected_time

        # SpinBox 업데이트 (값이 같으면 이벤트가 발생하지 않으므로 직접 업데이트)
        if self.time_spinbox.value() != selected_time:
            self.time_spinbox.setValue(selected_time)
        else:
            # 값이 같으면 수동으로 업데이트
            self.update_spectrum()
            self.update_time_line()
            self.update_texc_display()
            self.update_boltzmann_plot()  # 창3 업데이트

        # 마커 및 Texc 표시
        self.display_intensity_markers()

    def update_time_line(self):
        """창2의 시간 표시선 업데이트"""
        if self.time_line is not None:
            self.time_line.set_xdata([self.current_time, self.current_time])
            # 캔버스는 display_intensity_markers에서 그려지므로 여기서는 그리지 않음

    def calculate_boltzmann_plot_data(self):
        """
        Boltzmann Plot용 데이터 계산

        Returns:
        --------
        tuple : (energy_levels, y_values, slope, intercept, r_squared, texc_eV, error_msg)
            - energy_levels: X축 데이터 (E_n, eV)
            - y_values: Y축 데이터 (ln(I×λ/(g×A)))
            - slope: 선형 회귀 기울기
            - intercept: 선형 회귀 y절편
            - r_squared: 결정 계수
            - texc_eV: 여기 전자 온도 (eV)
            - error_msg: 오류 메시지 (성공 시 None)
        """
        if self.data is None or not self.balmer_timeseries:
            return None, None, None, None, None, None, "데이터가 로딩되지 않았습니다"

        # 현재 시간의 인덱스 찾기
        idx = self.get_time_index(self.current_time)

        energy_levels = []
        y_values = []
        intensities = {}

        for name in ['Hα', 'Hβ', 'Hγ']:
            # 인덱스 범위 검증
            if idx >= len(self.balmer_timeseries[name]):
                idx = len(self.balmer_timeseries[name]) - 1

            intensity = self.balmer_timeseries[name][idx]
            intensities[name] = intensity
            const = self.BALMER_CONSTANTS[name]

            # 유효하지 않은 intensity 체크
            if intensity <= 0:
                return None, None, None, None, None, None, "계산 불가: 유효하지 않은 Intensity"

            wavelength = const['wavelength']
            g = const['g']
            A = const['A']
            E = const['E']

            # Y값 계산: ln((I × λ) / (g × A))
            y = np.log((intensity * wavelength) / (g * A))

            energy_levels.append(E)
            y_values.append(y)

        # 선형 회귀
        energy_levels = np.array(energy_levels)
        y_values = np.array(y_values)

        slope, intercept, r_value, p_value, std_err = linregress(energy_levels, y_values)
        r_squared = r_value ** 2

        # 온도 계산: slope = -1 / (k_B × T)
        if slope >= 0:
            return None, None, None, None, None, None, "계산 불가: 비물리적 기울기 (양수)"

        T_kelvin = -1.0 / (self.K_B * slope)
        texc_eV = self.K_B * T_kelvin

        return energy_levels, y_values, slope, intercept, r_squared, texc_eV, None

    def on_boltzmann_checkbox_changed(self, state):
        """Boltzmann Plot 체크박스 상태 변경 이벤트"""
        if state == Qt.Checked:
            self.show_boltzmann_plot()
        else:
            self.close_boltzmann_plot()

    def show_boltzmann_plot(self):
        """Boltzmann Plot 창 표시"""
        if self.boltzmann_window is None:
            self.boltzmann_window = BoltzmannPlotWindow(self)

        # 데이터 계산 및 업데이트
        self.update_boltzmann_plot()

        # 창 표시
        self.boltzmann_window.show()
        self.boltzmann_window.raise_()
        self.boltzmann_window.activateWindow()

    def close_boltzmann_plot(self):
        """Boltzmann Plot 창 닫기"""
        if self.boltzmann_window is not None:
            self.boltzmann_window.close()
            self.boltzmann_window = None

    def update_boltzmann_plot(self):
        """Boltzmann Plot 창 업데이트"""
        if self.boltzmann_window is None or not self.boltzmann_window.isVisible():
            return

        # 데이터 계산
        energy_levels, y_values, slope, intercept, r_squared, texc_eV, error_msg = \
            self.calculate_boltzmann_plot_data()

        if error_msg is not None:
            # 계산 실패 시 메시지 표시
            self.boltzmann_window.show_error(error_msg)
        else:
            # 그래프 업데이트
            self.boltzmann_window.update_plot(
                energy_levels, y_values, slope, intercept,
                r_squared, texc_eV, self.current_time
            )


class BoltzmannPlotWindow(QWidget):
    """Boltzmann Plot 표시 팝업창 (Non-modal)"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_widget = parent
        self.setWindowTitle("Boltzmann Plot")
        self.setup_ui()

    def setup_ui(self):
        """UI 초기화"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        # Matplotlib Figure 생성
        self.figure = Figure(figsize=(6, 5), facecolor='white')
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)

        layout.addWidget(self.canvas)

        # 초기 빈 그래프 설정
        self.ax.set_xlabel('Upper Energy Level E$_n$ (eV)', fontsize=10)
        self.ax.set_ylabel('ln(I$_{nm}$λ$_{nm}$/g$_n$A$_{nm}$)', fontsize=10)
        self.ax.set_title('Boltzmann Plot', fontsize=10)
        self.ax.grid(True, linestyle='--', alpha=0.7, color='lightgray')
        self.figure.tight_layout()

        # 창 크기 설정
        self.resize(600, 500)

    def update_plot(self, energy_levels, y_values, slope, intercept, r_squared, texc_eV, current_time):
        """그래프 업데이트"""
        self.ax.clear()

        # 데이터 포인트 표시
        self.ax.plot(energy_levels, y_values, 'o',
                     markersize=10, color='#1f77b4', label='Data', zorder=3)

        # 선형 회귀 직선 표시
        x_min, x_max = energy_levels.min() - 0.3, energy_levels.max() + 0.3
        x_line = np.array([x_min, x_max])
        y_line = slope * x_line + intercept
        self.ax.plot(x_line, y_line, '--', color='#ff0000',
                     linewidth=1.5, label='Linear Fit', zorder=2)

        # 수치 표시 박스
        textstr = f'Texc = {texc_eV:.2f} eV\nSlope = {slope:.2f}\nR² = {r_squared:.2f}'
        self.ax.text(0.05, 0.05, textstr, transform=self.ax.transAxes,
                     fontsize=10, verticalalignment='bottom', horizontalalignment='left',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))

        # 축 라벨 및 제목
        self.ax.set_xlabel('Upper Energy Level E$_n$ (eV)', fontsize=10)
        self.ax.set_ylabel('ln(I$_{nm}$λ$_{nm}$/g$_n$A$_{nm}$)', fontsize=10)
        self.ax.set_title('Boltzmann Plot', fontsize=10)
        self.ax.grid(True, linestyle='--', alpha=0.7, color='lightgray')

        # 범례
        self.ax.legend(loc='upper right', fontsize=9)

        self.figure.tight_layout()
        self.canvas.draw()

        # 윈도우 제목 업데이트
        self.setWindowTitle(f'Boltzmann Plot - t = {current_time:.1f} s')

    def show_error(self, error_msg):
        """에러 메시지 표시"""
        self.ax.clear()
        self.ax.text(0.5, 0.5, error_msg, transform=self.ax.transAxes,
                     fontsize=12, verticalalignment='center', horizontalalignment='center',
                     bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='red'))
        self.ax.set_xlabel('Upper Energy Level E$_n$ (eV)', fontsize=10)
        self.ax.set_ylabel('ln(I$_{nm}$λ$_{nm}$/g$_n$A$_{nm}$)', fontsize=10)
        self.ax.set_title('Boltzmann Plot', fontsize=10)
        self.canvas.draw()

    def closeEvent(self, event):
        """창 닫힐 때 체크박스 동기화"""
        if self.parent_widget is not None:
            self.parent_widget.boltzmann_checkbox.setChecked(False)
        event.accept()


def main():
    """메인 함수"""
    app = QApplication(sys.argv)
    window = OESAnalyzer()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
