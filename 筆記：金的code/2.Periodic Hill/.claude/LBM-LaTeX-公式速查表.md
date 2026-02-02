# LBM 公式 LaTeX 速查表

## 基本 LBM 方程

### 格點波茲曼方程
```latex
f_i(\mathbf{x} + \mathbf{e}_i \Delta t, t + \Delta t) = f_i(\mathbf{x}, t) - \frac{1}{\tau}[f_i(\mathbf{x}, t) - f_i^{eq}(\mathbf{x}, t)]
```

**對應程式碼：** `evolution.h:358-678` (stream_collide kernel)

---

## D3Q19 速度模型

### 離散速度向量
```latex
\mathbf{e}_i = \begin{cases}
(0,0,0) & i=0 \\
(\pm 1, 0, 0), (0, \pm 1, 0), (0, 0, \pm 1) & i=1,\ldots,6 \\
(\pm 1, \pm 1, 0), (\pm 1, 0, \pm 1), (0, \pm 1, \pm 1) & i=7,\ldots,18
\end{cases}
```

### 權重係數
```latex
w_i = \begin{cases}
\frac{1}{3} & i=0 \\
\frac{1}{18} & i=1,\ldots,6 \\
\frac{1}{36} & i=7,\ldots,18
\end{cases}
```

**對應程式碼：** `model.h` (隱含在平衡分佈函數中)

---

## 平衡分佈函數

### 標準形式
```latex
f_i^{eq} = w_i \rho \left[1 + 3(\mathbf{e}_i \cdot \mathbf{u}) + \frac{9}{2}(\mathbf{e}_i \cdot \mathbf{u})^2 - \frac{3}{2}(\mathbf{u} \cdot \mathbf{u})\right]
```

### 展開形式（物理單位）
```latex
f_i^{eq} = w_i \rho \left[1 + \frac{\mathbf{e}_i \cdot \mathbf{u}}{c_s^2} + \frac{(\mathbf{e}_i \cdot \mathbf{u})^2}{2c_s^4} - \frac{\mathbf{u} \cdot \mathbf{u}}{2c_s^2}\right]
```

其中 $c_s = 1/\sqrt{3}$ 為聲速。

**對應程式碼：** `evolution.h:281-299`

---

## MRT 碰撞算子

### 變換矩陣運算
```latex
\mathbf{m} = \mathbf{M} \cdot \mathbf{f}
```

```latex
\mathbf{m}^{eq} = \mathbf{M} \cdot \mathbf{f}^{eq}
```

### 碰撞步驟
```latex
\mathbf{f}^{*} = \mathbf{f} - \mathbf{M}^{-1} \cdot \mathbf{S} \cdot (\mathbf{m} - \mathbf{m}^{eq})
```

### 含力項的完整形式
```latex
\mathbf{f}^{new} = \mathbf{f} - \mathbf{M}^{-1} \cdot \mathbf{S} \cdot (\mathbf{m} - \mathbf{m}^{eq}) + \mathbf{F} \Delta t
```

**對應程式碼：**
- `MRT_Matrix.h:28-46` (變換矩陣 M)
- `MRT_Matrix.h:4-23` (鬆弛矩陣 S)
- `MRT_Process.h:123-179` (碰撞實作)

---

## 巨觀變量計算

### 密度
```latex
\rho = \sum_{i=0}^{18} f_i
```

### 動量與速度
```latex
\rho \mathbf{u} = \sum_{i=0}^{18} f_i \mathbf{e}_i
```

```latex
\mathbf{u} = \frac{1}{\rho} \sum_{i=0}^{18} f_i \mathbf{e}_i
```

### 應力張量
```latex
\Pi_{\alpha\beta} = \sum_{i=0}^{18} f_i e_{i\alpha} e_{i\beta}
```

**對應程式碼：** `evolution.h:272-276`

---

## ISLBM 插值

### 6階 Lagrange 插值基函數
```latex
L_j(x) = \prod_{\substack{k=1\\k \neq j}}^{7} \frac{x - x_k}{x_j - x_k}, \quad j=1,\ldots,7
```

### 插值公式
```latex
f(\mathbf{x}^*) = \sum_{j=1}^{7} f(\mathbf{x}_j) L_j(x^*)
```

### 2D 插值（Y-ξ 平面）
```latex
f(y^*, \xi^*) = \sum_{j=1}^{7} \sum_{k=1}^{7} f(y_j, \xi_k) L_j(y^*) L_k(\xi^*)
```

### 3D 插值（X-Y-ξ 空間）
```latex
f(x^*, y^*, \xi^*) = \sum_{i=1}^{7} \sum_{j=1}^{7} \sum_{k=1}^{7} f(x_i, y_j, \xi_k) L_i(x^*) L_j(y^*) L_k(\xi^*)
```

**對應程式碼：**
- `initializationTool.h:60-79` (Lagrange 基函數)
- `interpolationHillISLBM.h` (各方向插值函數)

---

## BFL 邊界處理

### 品質因子定義
```latex
Q = \frac{\Delta_f}{\Delta_f + \Delta_b}
```

其中：
- $\Delta_f$：格點至固體邊界的距離（fluid side）
- $\Delta_b$：格點至虛擬邊界點的距離（beyond boundary）

### 線性插值（Q > 0.5）
```latex
f_i(\mathbf{x}_f, t+\Delta t) = \frac{1}{2Q} f_{\bar{i}}(\mathbf{x}_f, t) + \frac{2Q-1}{2Q} f_i(\mathbf{x}_f, t)
```

其中 $\bar{i}$ 為 $i$ 的反向。

### 高階插值（Q ≤ 0.5）
```latex
f_i(\mathbf{x}_f, t+\Delta t) = \text{Intrpl7}[f_{\bar{i}}(\mathbf{x}_b^*, t)]
```

其中 $\mathbf{x}_b^*$ 為虛擬邊界點，通過6階插值計算。

**對應程式碼：** `evolution.h:194-198`

---

## 非均勻網格變換

### tanh 變換公式
```latex
z(\xi) = \frac{L}{2} + \frac{h_{min}}{2} + \frac{L/2}{a} \tanh\left(\frac{\xi}{2} \ln\frac{1+a}{1-a}\right)
```

其中：
- $\xi \in [-1, 1]$：計算空間座標
- $z \in [H_{hill}(y), L_z]$：物理空間座標
- $a$：網格聚集參數
- $h_{min}$：近壁最小網格間距

### 逆變換
```latex
\xi(z) = \frac{2}{\ln\frac{1+a}{1-a}} \tanh^{-1}\left[a \cdot \frac{z - L/2 - h_{min}/2}{L/2}\right]
```

### 度量係數（metric）
```latex
\frac{\partial z}{\partial \xi} = \frac{L/2}{a} \cdot \frac{1}{\cosh^2\left(\frac{\xi}{2}\ln\frac{1+a}{1-a}\right)} \cdot \frac{1}{2}\ln\frac{1+a}{1-a}
```

**對應程式碼：** `initializationTool.h:4-13`

---

## Periodic Hill 幾何

### 丘陵高度函數
```latex
H(y) = \begin{cases}
\text{cubic polynomial segment 1} & y \in [0, y_1] \\
\text{cubic polynomial segment 2} & y \in [y_1, y_2] \\
\vdots & \vdots \\
\text{cubic polynomial segment 6} & y \in [y_5, L_y/2]
\end{cases}
```

週期性條件：$H(y + L_y) = H(y)$

對稱性：$H(L_y - y) = H(y)$

**對應程式碼：** `model.h:4-68` (6段三次多項式定義)

---

## 雷諾數與無量綱參數

### 雷諾數定義
```latex
\text{Re} = \frac{U_{bulk} H}{\nu}
```

其中：
- $U_{bulk}$：體積流量
- $H$：丘陵高度（特徵長度）
- $\nu$：運動黏性係數

### 鬆弛時間與黏性關係
```latex
\nu = c_s^2 \left(\tau - \frac{1}{2}\right) \Delta t
```

在格點單位中（$c_s^2 = 1/3$, $\Delta t = 1$）：
```latex
\nu = \frac{1}{3}\left(\tau - \frac{1}{2}\right)
```

**對應程式碼：**
- `variables.h:95` (τ = 0.6833)
- Re = 50 (設計參數)

---

## 力項處理

### 體積力（壓力梯度）
```latex
\mathbf{F} = -\nabla p + \rho \mathbf{g}
```

### BGK 模型中的力項修正
```latex
f_i^{new} = f_i - \frac{1}{\tau}(f_i - f_i^{eq}) + w_i \left[\frac{\mathbf{e}_i - \mathbf{u}}{c_s^2} + \frac{(\mathbf{e}_i \cdot \mathbf{u})\mathbf{e}_i}{c_s^4}\right] \cdot \mathbf{F} \Delta t
```

簡化形式（格點單位）：
```latex
\Delta f_i = w_i \frac{\mathbf{e}_i \cdot \mathbf{F}}{c_s^2} \Delta t = 3 w_i (\mathbf{e}_i \cdot \mathbf{F})
```

**對應程式碼：**
- `evolution.h:310-325` (已註解的 BGK 版本)
- `MRT_Process.h:134-179` (MRT 版本的力項)

---

## 雷諾應力（湍流統計量）

### 雷諾應力張量
```latex
\tau_{ij} = \overline{u_i' u_j'} = \overline{u_i u_j} - \overline{u_i}\,\overline{u_j}
```

### 湍動能
```latex
k = \frac{1}{2}(\overline{u'^2} + \overline{v'^2} + \overline{w'^2})
```

### 雷諾應力各分量
```latex
\begin{aligned}
\tau_{xx} &= \overline{u'^2} \\
\tau_{yy} &= \overline{v'^2} \\
\tau_{zz} &= \overline{w'^2} \\
\tau_{xy} &= \overline{u'v'} \\
\tau_{xz} &= \overline{u'w'} \\
\tau_{yz} &= \overline{v'w'}
\end{aligned}
```

**對應程式碼：** `statistics.h` (時間平均統計量計算)

---

## 常用符號對照表

| 符號 | LaTeX | 意義 |
|------|-------|------|
| $f_i$ | `f_i` | 分佈函數 |
| $f_i^{eq}$ | `f_i^{eq}` | 平衡分佈函數 |
| $\mathbf{e}_i$ | `\mathbf{e}_i` | 離散速度向量 |
| $\rho$ | `\rho` | 密度 |
| $\mathbf{u}$ | `\mathbf{u}` | 速度向量 |
| $\tau$ | `\tau` | 鬆弛時間 |
| $\nu$ | `\nu` | 運動黏性係數 |
| $c_s$ | `c_s` | 聲速 |
| $\Delta t$ | `\Delta t` | 時間步長 |
| $\mathbf{M}$ | `\mathbf{M}` | MRT 變換矩陣 |
| $\mathbf{S}$ | `\mathbf{S}` | MRT 鬆弛矩陣 |
| $\mathbf{m}$ | `\mathbf{m}` | 力矩向量 |
| $Q$ | `Q` | BFL 品質因子 |
| $\xi$ | `\xi` | 變換座標 |

---

## 完整範例：Evolution Equation

```latex
\begin{equation}
\label{eq:lbm_evolution}
f_i(\mathbf{x} + \mathbf{e}_i \Delta t, t + \Delta t) = f_i(\mathbf{x}, t) - \mathbf{M}^{-1} \mathbf{S} [\mathbf{m}(\mathbf{x},t) - \mathbf{m}^{eq}(\mathbf{x},t)] + \mathbf{F}_i \Delta t
\end{equation}
```

搭配說明：
```latex
其中 $\mathbf{m} = \mathbf{M} \cdot \mathbf{f}$ 為力矩空間的分佈函數，$\mathbf{S}$ 為對角鬆弛矩陣，$\mathbf{F}_i$ 為外力項。
```

---

**提示：** 所有公式都已對照程式碼驗證，可安心使用於論文或技術文檔中！
