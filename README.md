# OES 데이터 분석 프로그램

수소 플라즈마 OES(Optical Emission Spectroscopy) 데이터를 분석하고 시각화하는 Python GUI 프로그램입니다.

## 기능

- **스펙트럼 분석**: 특정 시간에서의 파장별 발광 강도 시각화
- **시계열 분석**: 발머 계열 파장(Hα, Hβ, Hγ)의 시간에 따른 강도 변화 추적
- **가우시안 가중 평균**: 정확한 파장 강도 추출을 위한 고급 알고리즘
- **여기 전자 온도 계산**: Boltzmann Plot Method 기반 Texc 계산 및 표시
- **Boltzmann Plot 시각화**: 온도 계산 과정을 그래프로 시각화하는 팝업창
- **인터랙티브 UI**: 그래프 클릭 시 Intensity 마커 및 Texc 실시간 표시

## 요구사항

```bash
pip install numpy pandas matplotlib PyQt5 scipy
```

## 실행 방법

```bash
python oes_analysis.py
```

## 사용 방법

1. **파일 로딩**
   - "파일 로딩" 버튼을 클릭하여 .dat 파일 선택
   - 샘플 데이터: `sample_oes_data.dat`

2. **Window 설정**
   - 가우시안 윈도우 크기 조정 (기본값: 2.0 nm)
   - 범위: 0.5 ~ 10.0 nm

3. **시간 선택**
   - SpinBox를 사용하여 직접 입력
   - 또는 하단 그래프(시계열)를 클릭하여 선택
   - 클릭 시 해당 시간의 Intensity 마커 및 Texc 값 표시

4. **그래프 해석**
   - **상단 그래프**: 선택된 시간의 전체 스펙트럼 (200-800 nm)
   - **하단 그래프**: 발머 계열 파장의 시간에 따른 강도 변화 (좌측 Y축) 및 Texc 변화 (우측 Y축)
     - Hα (656.28 nm): 수소 알파선 (파란색)
     - Hβ (486.13 nm): 수소 베타선 (주황색)
     - Hγ (434.05 nm): 수소 감마선 (초록색)
     - Texc: 여기 전자 온도 (검정색, 사각 마커)
   - 빨간 점선: 현재 선택된 시간 위치
   - 클릭 시 각 라인의 Intensity 값 및 Texc 값 표시

5. **여기 온도 확인**
   - 컨트롤 패널: 현재 시간의 Texc 및 R² 값
   - 그래프 우측 상단: 클릭 시 해당 시간의 Texc 및 R² 값

6. **Boltzmann Plot 시각화**
   - 컨트롤 패널의 "Boltzmann Plot" 체크박스를 선택하여 팝업창 활성화
   - 팝업창은 비모달(non-modal) 방식으로 동작하여 메인 창과 동시에 조작 가능
   - 시간 변경 시 자동으로 그래프 업데이트
   - 그래프 구성:
     - 3개의 파란색 원형 마커: Hα, Hβ, Hγ 데이터 포인트
     - 빨간색 점선: 선형 회귀 라인 (Linear Fit)
     - X축: Upper Energy Level (eV)
     - Y축: ln(I × λ / (g × A))
   - 텍스트 박스: Texc (eV), Slope, R² 값 표시
   - 창 제목: 현재 시간 표시 (예: "Boltzmann Plot - t = 1.5 s")

## 입력 파일 형식

- 확장자: `.dat` (탭으로 구분된 텍스트 파일)
- 열 1: 센서명 (문자열)
- 열 2: Run Time (초 단위, float)
- 열 3-1203: 파장별 광강도 (200.0-800.0 nm, 0.5 nm 간격)

예시:
```
Name    Run Time    200.0    200.5    201.0    ...    800.0
A       0.5         1234     1235     1240     ...    890
A       1.0         1250     1260     1255     ...    900
```

## 핵심 알고리즘

### 가우시안 가중 평균

각 발머 계열 파장의 강도는 다음과 같이 계산됩니다:

```
w(λ) = exp(-((λ - λ₀)² / (2σ²)))
σ = Window / 2
Intensity = Σ(w(λ) × I(λ)) / Σ(w(λ))
```

여기서:
- `λ₀`: 중심 파장 (Hα, Hβ, Hγ)
- `Window`: 사용자가 설정한 윈도우 크기
- `I(λ)`: 파장 λ에서의 측정 강도

### Boltzmann Plot Method (여기 전자 온도 계산)

여기 전자 온도(Texc)는 다음과 같이 계산됩니다:

```
ln(I × λ / (g × A)) = -E / (k_B × T_exc) + C
```

계산 절차:
1. 각 발머 라인(Hα, Hβ, Hγ)에 대해 y = ln(I × λ / (g × A)) 계산
2. x = E_n (상위 에너지 레벨) 설정
3. 선형 회귀로 기울기(slope) 계산
4. T_exc = -1 / (k_B × slope) 로 온도(K) 계산
5. eV 단위로 변환: T_exc_eV = k_B × T_exc

여기서:
- `I`: 측정된 발광 강도 (가우시안 가중 평균)
- `λ`: 파장 (nm)
- `g`: 상위 에너지 준위의 통계적 가중치 (degeneracy)
- `A`: 전이 확률 (s⁻¹)
- `E`: 상위 에너지 레벨 (eV)
- `k_B`: 볼츠만 상수 (8.617333262 × 10⁻⁵ eV/K)

NIST 원자 스펙트럼 상수 (하드코딩):
- Hα (656.28 nm): g=18, A=4.41×10⁷ s⁻¹, E=12.09 eV
- Hβ (486.13 nm): g=32, A=8.42×10⁶ s⁻¹, E=12.75 eV
- Hγ (434.05 nm): g=50, A=2.53×10⁶ s⁻¹, E=13.05 eV

## 프로그램 구조

```
oes_analysis.py
├── OESAnalyzer (메인 클래스)
│   ├── init_ui(): GUI 초기화
│   ├── create_control_panel(): 컨트롤 패널 생성 (온도 표시 및 Boltzmann Plot 체크박스 포함)
│   ├── create_graph_area(): 그래프 영역 생성 (툴바 제거)
│   ├── load_file(): 데이터 파일 로딩
│   ├── gaussian_weighted_average(): 가우시안 가중 평균 계산
│   ├── calculate_texc(): Boltzmann Plot Method로 Texc 계산
│   ├── calculate_boltzmann_plot_data(): Boltzmann Plot 데이터 계산
│   ├── update_spectrum(): 스펙트럼 그래프 업데이트
│   ├── update_timeseries(): 시계열 그래프 업데이트 (보조 Y축 포함)
│   ├── update_texc_display(): 컨트롤 패널의 Texc 표시 업데이트
│   ├── display_intensity_markers(): 창2에 Intensity 마커 및 Texc 표시
│   ├── clear_intensity_markers(): 마커 및 annotation 제거
│   ├── on_timeseries_click(): 창2 클릭 이벤트 처리
│   ├── on_boltzmann_checkbox_changed(): Boltzmann Plot 체크박스 상태 변경 처리
│   ├── show_boltzmann_plot(): Boltzmann Plot 창 표시
│   ├── close_boltzmann_plot(): Boltzmann Plot 창 닫기
│   └── update_boltzmann_plot(): Boltzmann Plot 창 업데이트
├── BoltzmannPlotWindow (팝업 창 클래스)
│   ├── __init__(): 창 초기화 및 레이아웃 생성
│   ├── update_plot(): 그래프 및 텍스트 업데이트
│   └── closeEvent(): 창 닫기 이벤트 처리 (체크박스 동기화)
└── main(): 메인 함수
```

## GUI 구성

### 컨트롤 패널 (좌측)
- 파일 로딩 버튼
- Window 입력 (0.5-10.0 nm)
- 시간 입력 (데이터 범위)
- **여기 온도 표시** (읽기 전용, 복사 가능)
- **R² 표시**
- **Boltzmann Plot 체크박스**: 시각화 팝업창 표시/숨김

### 그래프 영역 (우측)
- **창1 (상단)**: 스펙트럼 그래프
  - X축: Wavelength (nm)
  - Y축: Emission Intensity (a.u.)

- **창2 (하단)**: 시계열 그래프
  - 좌측 Y축: Intensity (a.u.) - Hα, Hβ, Hγ
  - 우측 Y축: Excitation Temperature (eV) - Texc
  - 클릭 시 마커 및 Texc 값 표시

- **창3 (팝업)**: Boltzmann Plot 그래프
  - X축: Upper Energy Level (eV)
  - Y축: ln(I × λ / (g × A))
  - 데이터 포인트: 파란색 원형 마커
  - 선형 회귀: 빨간색 점선
  - 정보 표시: Texc, Slope, R² 값

## 개발자 정보

- Python 3.6 이상 필요
- PEP8 스타일 가이드 준수
- 한글 주석 포함

## 라이선스

MIT License
