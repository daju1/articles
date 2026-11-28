# MenDrive_qnm.ipynb: объёмная сила на стенку против тензора Максвелла — проверка (2026-11-28)

Объект: `rezonatory_i_volnovody/MenDrive_qnm.ipynb` (909 ячеек; `driving_mode = driving_mode_1 = 'je_y, k_y=0, TE'`, ток $j_{e,y}^{ext}$ в левой стенке),
библиотека `rezonatory_i_volnovody/mendrivepylib/` и нелинейная 1D-модель `mendrive_llg/MenDriveCpp/standalone/`.
Исходные файлы не изменялись. Числа ниже получены независимым пересчётом на `numpy` (`longdouble`, 33 знака) в облаке.
Тот же расчёт оформлен как проверочный ноутбук `check_force_tensor_vs_volume.ipynb`: 17 ячеек, все проверки OK.

Обозначения: СГС, $e^{-i\omega t}$, $e^{ik_zz}$, $k_z=\pi/L_z$, $L_z=100$ см, $k_y=0$. Левая стенка (металл) занимает $x<-a$, вакуум $|x|\le a$, правая стенка (феррит) $x>a$, $a=0.5$ см.
$\varepsilon'=\varepsilon+4\pi i\sigma_e/\omega$, $\mu'=\mu+4\pi i\sigma_m/\omega$. Усреднение по периоду: $\langle AB\rangle=\tfrac12\mathrm{Re}(AB^*)$.
Параметры (вывод ячейки 420): левая стенка $\varepsilon=\mu=1$, $\sigma_e=5.625\cdot10^{17}\ \text{с}^{-1}$, $\sigma_m=0$. Правая: $\varepsilon=16$, $\mu=100$, $\sigma_e=9.9\cdot10^{12}$, $\sigma_m=6.48\cdot10^{7}\ \text{с}^{-1}$.
$h_l=1.875\cdot10^{-4}$ см, $h_r=6.25\cdot10^{-3}$ см, $j_0=10^6$. QNM: $\omega=1.79881421\cdot10^{11}-7.5654849\cdot10^{9}i$.

---

## 0. ТЕЗИС

1. **3466 Н/кВт — это один член, а не сила.** На 99.98 % это $F_{I\times\mathrm{rot}B}$ (ячейки 293/306/378), посчитанный с **потерянным комплексным сопряжением** $\overline{\mu'_r}\to\mu'_r$.
   Символы `mu_r_xx_c, mu_r_yy_c, mu_r_zz_c` в `mendrivepylib/variables/common.py` не объявлены `complex` (объявлены только `mu_l_*_c` и `epsilon_*_c`).
   После `.simplify()` в ячейках 291 и 293 Maxima считает их вещественными и выбрасывает `conjugate`: в выводе ячеек 291/293/303 справа стоит `mu_r_zz_c^2` и `mu_r_xx_c` без `conjugate`, слева `conjugate(mu_l_zz_c)`.
   Поскольку $\mathrm{Im}\,\mu'_r=4\pi\sigma_m/\omega=814$ при $\omega=10^6$, вместо $|\mu'|^2$ получается $\mu'^2$:
   $$F^{bug}_{I\times\mathrm{rot}B}\simeq+\frac{|H_z(a)|^2}{16\pi}\Big(\frac{4\pi\sigma_m}{\omega}\Big)^2 \;\propto\;\frac{1}{\omega^2}.$$
   Это и есть зависимость $1/\omega^2$: она идёт от магнитной проводимости $\sigma_m$, а не от резонатора.
2. Даже с верным сопряжением этот член и ещё 8 из 12 слагаемых `force_dyn_right` (ячейка 396) — **двойной счёт**: «Абрагамовы» $P\times\mathrm{rot}E$, $I\times\mathrm{rot}B$, динамические $\dot P\times B$, $\dot M\times D$ и «смайтовские» поверхностные члены.
   Кроме того, $\mathbf P'=\frac{\varepsilon'-1}{4\pi}\mathbf E$ и $\mathbf M'=\frac{\mu'-1}{4\pi}\mathbf H$ содержат токи проводимости. Правильная сила на правую стенку:
   $F_r=F_{jmD}+p_{Bt,jeB}+[T_{xx}]_a=-\langle T^{vac}_{xx}\rangle(a)$. Это и есть `N_per_kW_maxwell_right_wall` $=-0.0071$ (воспроизведено до $10^{-15}$).
3. Поля в `scan_omega_re` — **не решение**. `calc_thrust` берёт минимальный сингулярный вектор однородной $M_4$ (`svd_X4`) и не использует $B_{ext}$.
   `residual≈1e-3` мал только потому, что строки сшивки $E_y$ масштабированы как $\omega/(ck_z)\sim10^{-3}$; относительная невязка сшивки $E_y$ составляет 0.3–1.6.
4. В вынужденной задаче (ячейки 21/75/232/477) источник с профилем собственной моды требует **секулярного члена**. В граничном условии должен стоять $2\pi$, а не $4\pi$: все амплитуды ячейки 477 завышены ровно в 2 раза.
5. При правильном расчёте **объём = тензор** до $10^{-16}$ (Минковский и Эйнштейн–Лауб), $P_{src}=P_{diss}$ до $10^{-18}$, а **суммарная сила равна нулю** ($|F_l+F_r|/|F_r|\sim10^{-11}$).
   Стенки давят друг на друга с силой $\sim10^{-7}$ Н/кВт при $\omega=10^6\ldots10^7$ и $2\cdot10^{-9}$ Н/кВт на частоте QNM.

---

## 1. Что считает ноутбук (формулы)

### 1.1 Поля (ячейки 45–48, 101–131; $k_y=0$, `sign_K=1`, $\mathrm{Im}K>0$ — `fix_K_sign`)

$$\text{вакуум: } H_z=A_1e^{-iK_vx}+A_2e^{iK_vx},\quad H_x=\frac{c^2k_zK_v}{c^2k_z^2-\omega^2}\big(A_2e^{iK_vx}-A_1e^{-iK_vx}\big),\quad E_y=\frac{cK_v\omega}{c^2k_z^2-\omega^2}\big(A_1e^{-iK_vx}-A_2e^{iK_vx}\big),$$
$$\text{левая: } H_z=B_{H,zl}e^{-iK_l(x+a)},\quad E_y=\frac{cK_l\mu'_l\omega B_{H,zl}+4\pi i\mu'_l\omega j_0}{c^2k_z^2-\varepsilon'_l\mu'_l\omega^2}e^{-iK_l(x+a)},\qquad
\text{правая: } H_z=B_{H,zr}e^{iK_r(x-a)},\quad E_y=-\frac{cK_r\mu'_r\omega B_{H,zr}}{c^2k_z^2-\varepsilon'_r\mu'_r\omega^2}e^{iK_r(x-a)},$$
$$K_v=\sqrt{\omega^2/c^2-k_z^2},\qquad K^2=\Big(\varepsilon'\mu'\frac{\omega^2}{c^2}-k_z^2\Big)\frac{\mu'_{zz}}{\mu'_{xx}}\ \ (\text{ячейки 171–176}),\qquad \mathbf B=\mu'\mathbf H,\ \mathbf D=\varepsilon'\mathbf E\ \ (\text{ячейки 51, 70}).$$
Сшивка $H_z$ и $E_y$ при $x=\mp a$ даёт $M_4X=B_{ext}j_0$, $X=(B_{H,zl},A_1,A_2,B_{H,zr})$, где $B_{ext}[2]=-\frac{4\pi i\mu'_l\omega}{c^2k_z^2-\varepsilon'_l\mu'_l\omega^2}$ (ячейки 222–232, 473).

### 1.2 Объёмные силы (ячейки 291–347), $\int_a^\infty$ / $\int_{-\infty}^{-a}$, $\int e^{-2\mathrm{Im}K|x\mp a|}=\frac{1}{2\mathrm{Im}K}$

$$F_{jmD}=\frac{1}{2c}\mathrm{Re}\!\int\sigma_mH_z(\varepsilon'E_y)^*,\qquad p_{Bt,jeB}=\frac{1}{2c}\mathrm{Re}\!\int\sigma_eE_y(\mu'H_z)^*,\qquad p_{Ht,jeH}=\frac{1}{2c}\mathrm{Re}\!\int\sigma_eE_yH_z^*\ (\text{не входит в сумму}),$$
$$F_{I\times\mathrm{rot}B}=-\frac12\mathrm{Re}\!\int\frac{\mu'-1}{4\pi}H_z\frac{i\omega}{c}\big(\mu'\varepsilon'E_y\big)^*,\qquad
F_{dyn,I,D}=\frac12\mathrm{Re}\!\int\frac{i\omega}{c}\frac{\mu'-1}{4\pi}H_z(\varepsilon'E_y)^*,$$
$$F_{P\times\mathrm{rot}E}=-\frac12\mathrm{Re}\!\int\frac{\varepsilon'-1}{4\pi}E_y\frac{i\omega}{c}(\mu'H_z)^*,\qquad
F_{dyn,P,B}=\frac12\mathrm{Re}\!\int\Big(-\frac{i\omega}{c}\Big)\frac{\varepsilon'-1}{4\pi}E_y(\mu'H_z)^*,\qquad
p_{Bt,je\,ext}=\frac{1}{2c}\mathrm{Re}\!\int j^{ext}_y(\mu'_lH_z)^*.$$

### 1.3 Поверхностные «смайтовские» члены (ячейки 381–390)

$$p_{I,conv,Ht}=\frac{\mathrm{Re}(1-\mu'_{zz})|H^{vac}_z(a)|^2}{16\pi},\quad t_{I,conv,Hn}=\frac{\mathrm{Re}(1-\mu'_{xx})|H_x(a^+)|^2}{16\pi},\quad
p_{Et}=\frac{\mathrm{Re}(1-\varepsilon'_{yy})|E^{vac}_y(a)|^2}{16\pi},\quad t_{Dn}=\frac{\mathrm{Re}(1-\varepsilon'_{xx})|E_x(a^+)|^2}{16\pi}=0.$$

### 1.4 Тензор (`mendrivepylib/forces/common.py`, `compute_maxwell_stress_tensor_symbolic`, Тамм 105.1)

$$\langle T_{xx}\rangle=\frac{1}{16\pi}\mathrm{Re}\big[2E_xD_x^*+2H_xB_x^*-\mathbf E\cdot\mathbf D^*-\mathbf H\cdot\mathbf B^*\big]$$
(при $\mathbf D=\varepsilon'\mathbf E$ мнимая часть $\varepsilon'$ под $\mathrm{Re}$ выпадает). Ячейка 396:
`f_jump_maxwell_right` $=T_r(a)-T_v(a)$, `f_maxwell_right_wall` $=-T_r(a)+[T_r(a)-T_v(a)]=-T_v(a)$, `f_maxwell_left_wall` $=T_v(-a)$.

### 1.5 Сумма и нормировка (ячейка 396)

$$F^{nb}_{right}=F_{jmD}+p_{Bt,jeB}+p_{I,conv,Ht}+t_{I,conv,Hn}+F_{I\times\mathrm{rot}B}+F_{dyn,I,D}+t_{Dn}+p_{Et}+F_{P\times\mathrm{rot}E}+F_{dyn,P,B}+[T_{xx}]_a+F_{jm,ext,D},$$
$$W=\int_{-a-h_l}^{-a}\!+\!\int_a^{a+h_r}\tfrac12\big(\sigma_e|\mathbf E|^2+\sigma_m|\mathbf H|^2\big)dx\quad(\text{ячейка 392; только слои }h),\qquad
\text{Н/кВт}=\frac{F[\text{дин/см}^2]\cdot10^{-5}}{W[\text{эрг/(с·см}^2)]\cdot10^{-7}}\cdot10^{3}=10^5\frac{F}{W}.$$
Мощность источника $P_{src}$ не вычисляется. Поток Пойнтинга тоже не используется.

---

## 2. Правильный баланс для линейной среды

Свободные токи $\mathbf j_e=\sigma_e\mathbf E+\mathbf j^{ext}$, $\mathbf j_m=\sigma_m\mathbf H$; $\varepsilon,\mu$ вещественные и кусочно-постоянные. Уравнения:
$\mathrm{rot}\,\mathbf H=\frac1c\dot{\mathbf D}+\frac{4\pi}{c}\mathbf j_e$, $\mathrm{rot}\,\mathbf E=-\frac1c\dot{\mathbf B}-\frac{4\pi}{c}\mathbf j_m$.

**Минковский.** Внутри однородной области $\partial_jT^M_{ij}-\partial_tg^M_i=f^M_i$, где
$f^M=\rho_e\mathbf E+\rho_m\mathbf H+\frac1c\mathbf j_e\times\mathbf B-\frac1c\mathbf j_m\times\mathbf D$.
Для TE с $k_y=0$ $\nabla\!\cdot\mathbf E=0$, а $\nabla\!\cdot(\mu'\mathbf H)=0\Rightarrow\nabla\!\cdot\mathbf H=0$, поэтому $\rho_e=\rho_m=0$. В стационарном периодическом режиме ($\omega$ вещественная) $\langle\partial_t\mathbf g\rangle=0$, и
$$\partial_x\langle T^M_{xx}\rangle=\langle f^M_x\rangle=\frac{1}{2c}\mathrm{Re}\big(j_{e,y}B_z^*\big)+\frac{1}{2c}\mathrm{Re}\big(j_{m,z}D_y^*\big),\qquad
\int_a^\infty\langle f^M_x\rangle dx=-\langle T^M_{xx}\rangle(a^+),\quad \int_{-\infty}^{-a}\langle f^M_x\rangle dx=\langle T^M_{xx}\rangle(-a^-).$$
Поверхностная сила (скачок $\varepsilon,\mu$) равна $[T^M]$. Полная сила на тело, граничащее с вакуумом:
$$F_r=-\langle T^{vac}_{xx}\rangle(a),\qquad F_l=\langle T^{vac}_{xx}\rangle(-a),\qquad F_l+F_r=-\int_{-a}^a\partial_x\langle T^{vac}_{xx}\rangle dx=0,$$
так как в вакууме $\partial_x\langle T_{xx}\rangle=-\langle\partial_tg_x\rangle=0$. Этот результат **не зависит от формулировки** (Минковский / Абрагам / Эйнштейн–Лауб): в вакууме все тензоры совпадают.
Он верен и для нелинейных стенок в любом периодическом режиме.

*Замечание о $\varepsilon',\mu'$ в ноутбуке.* В $F_{jmD}$ стоит $D'=\varepsilon'E$, в $p_{Bt,jeB}$ стоит $B'=\mu'H$. Каждый содержит паразитный член $\mp\frac{4\pi\sigma_e\sigma_m}{\omega}\mathrm{Re}(iE_yH_z^*)/2c$, но в сумме они **сокращаются**: проверено, $F_{jmD}+p_{Bt,jeB}=-T^M(a^+)$ с точностью $10^{-16}$.

**Эйнштейн–Лауб** (независимая проверка объёмной плотности, $\mathbf P=\frac{\varepsilon-1}{4\pi}\mathbf E$, $\mathbf M=\frac{\mu-1}{4\pi}\mathbf H$):
$$\langle f^{EL}_x\rangle=\langle(\mathbf M\!\cdot\!\nabla)H_x\rangle+\frac{1}{2c}\mathrm{Re}\big[(j_e-i\omega P)_yH_z^*\big]+\frac{1}{2c}\mathrm{Re}\big[(j_m-i\omega M)_zE_y^*\big],\qquad
\langle T^{EL}_{xx}\rangle=\frac{1}{16\pi}\mathrm{Re}\big[2B_xH_x^*-|\mathbf E|^2-|\mathbf H|^2\big].$$

**Двойной счёт в ноутбуке.**
$\langle\dot{\mathbf P}\times\mathbf B\rangle=-\langle\mathbf P\times\dot{\mathbf B}\rangle=c\langle\mathbf P\times\mathrm{rot}\,\mathbf E\rangle$, следовательно $F_{P\times\mathrm{rot}E}\equiv F_{dyn,P,B}$ (проверено: равны до $10^{-16}$).
$\dot{\mathbf P}'=-i\omega\frac{\varepsilon-1}{4\pi}\mathbf E+\sigma_e\mathbf E$, то есть ток проводимости посчитан ещё раз. В левой стенке ($\varepsilon=1$) $F_{P\times\mathrm{rot}E,l}=F_{dyn,P,B,l}=p_{Bt,jeB,l}$: тройной счёт.
Аналогично $\mathbf M'=\frac{\mu'-1}{4\pi}\mathbf H=\frac{\mu-1}{4\pi}\mathbf H+\frac{i\sigma_m}{\omega}\mathbf H$ содержит магнитный ток. Отсюда $F^{ok}_{I\times\mathrm{rot}B}\simeq-\frac{|H_z|^2}{16\pi}\big[\mu(\mu-1)+\mu''^2\big]$.
Поверхностные $p_{I,conv,Ht},t_{I,conv,Hn},p_{Et},t_{Dn}$ — части $[T^M]_a$, а $[T^M]_a$ уже включён как `f_jump_maxwell_right`.

---

## 3. Численная проверка

### 3.1 Воспроизведение ячейки 479 и источник 3466 / $1/\omega^2$

Поля `svd_X4`, формулы ячеек 291–396. `thrust_N_per_kW_right` воспроизведено до $10^{-14}$, `thrust_N_per_kW_left` — до $10^{-14}$ (с учётом ошибки E6), `residual` — до $10^{-16}$.

| $\omega$, рад/с | ноутбук: правая стенка, Н/кВт | ноутбук: $-T^{vac}_{xx}(a)$, Н/кВт | воспроизведение (12 членов) | из них $F_{I\times\mathrm{rot}B}$ (баг $\overline{\mu'}\to\mu'$) | тот же член с верным $\overline{\mu'}$ | сумма 12 членов с верным $\overline{\mu'}$ | $(4\pi\sigma_m/\omega)^2$ | `residual` SVD | $W_h/W_{полн}$ |
|---|---|---|---|---|---|---|---|---|---|
| 1.0000e+06 | 3466.33 | -0.007122 | 3466.33 | 3465.55 | -3516.9 | -3516.65 | 6.631e+05 | 0.00098 | 0.124 |
| 2.0000e+06 | 857.47 | 0.002949 | 857.47 | 856.688 | -908.104 | -907.831 | 1.658e+05 | 0.00112 | 0.125 |
| 3.0000e+06 | 373.963 | 0.004158 | 373.963 | 373.188 | -424.617 | -424.335 | 7.368e+04 | 0.00117 | 0.126 |
| 4.0000e+06 | 204.351 | 0.004457 | 204.351 | 203.587 | -254.982 | -254.688 | 4.144e+04 | 0.00121 | 0.128 |
| 5.0000e+06 | 125.512 | 0.00451 | 125.512 | 124.762 | -176.051 | -175.746 | 2.652e+04 | 0.00125 | 0.129 |
| 6.0000e+06 | 82.4241 | 0.004465 | 82.4241 | 81.6898 | -132.783 | -132.468 | 1.842e+04 | 0.0013 | 0.131 |
| 7.0000e+06 | 56.2515 | 0.004375 | 56.2515 | 55.5339 | -106.337 | -106.011 | 1.353e+04 | 0.00135 | 0.133 |
| 8.0000e+06 | 39.1341 | 0.004261 | 39.1341 | 38.4339 | -88.8505 | -88.5163 | 1.036e+04 | 0.0014 | 0.135 |
| 9.0000e+06 | 27.3175 | 0.004137 | 27.3175 | 26.635 | -76.5776 | -76.2359 | 8186 | 0.00145 | 0.138 |
| 1.0000e+07 | 18.8217 | 0.004008 | 18.8217 | 18.157 | -67.5486 | -67.201 | 6631 | 0.00151 | 0.14 |
| 1.7988e+11 | — | — | -0.0594783 | -0.0606828 | -0.0606884 | -0.059484 | 2.049e-05 | 0.114 | 1 |

Асимптотика при $\mathrm{rot}H_y\simeq-iK_rH_z$, $K_r=K'+iK''$, $\mu''=4\pi\sigma_m/\omega$ (совпадает с точным членом до $10^{-5}$):
$$F^{bug}=-\frac{|H_z(a)|^2}{16\pi}\Big[\mu(\mu-1)-\mu''^2-\mu''(2\mu-1)\frac{K'}{K''}\Big],\qquad F^{ok}=-\frac{|H_z(a)|^2}{16\pi}\Big[\mu(\mu-1)+\mu''^2-\mu''\frac{K'}{K''}\Big].$$

Разбор `force_dyn_right` при $\omega=10^6$:

| слагаемое `force_dyn_right` | Н/кВт (как в ноутбуке) | в правильной формулировке |
|---|---|---|
| $F_{jmD}$ (яч. 320–325, 361) | 0.260335 | входит |
| $p_{Bt,jeB}$ (яч. 339–347, 366) | 0.262286 | входит |
| $p_{I,conv,Ht}=\frac{(1-\mu)|H^{vac}_z|^2}{16\pi}$ (яч. 382) | -0.517401 | лишний (двойной счёт) |
| $t_{I,conv,Hn}$ (яч. 384) | -4.49658e-06 | лишний (двойной счёт) |
| $F_{I\times\mathrm{rot}B}$ (яч. 293, 306, 378) | 3465.55 | лишний (двойной счёт) |
| $F_{dyn,I,D}$ (яч. 314–318, 356) | 0.257059 | лишний (двойной счёт) |
| $t_{Dn}$ (яч. 390) | 0 | лишний (двойной счёт) |
| $p_{Et}$ (яч. 388) | -2.0881e-07 | лишний (двойной счёт) |
| $F_{P\times\mathrm{rot}E}$ (яч. 291, 300, 379) | 0.782956 | лишний (двойной счёт) |
| $F_{dyn,P,B}$ (яч. 308–312, 358) | 0.262286 | лишний (двойной счёт) |
| $[T_{xx}]_a$ `f_jump_maxwell_right` | -0.529743 | входит |
| $F_{jm,ext,D}$ (яч. 327–332) | 0 | входит |
| **сумма** | **3466.33** | $F_{jmD}+p_{Bt,jeB}+[T]_a=-0.00712183=-T^{vac}_{xx}(a)$ |

### 3.2 Правильный расчёт (вынужденное решение, секулярный источник, $j_0=10^6$, полубесконечные стенки, $k_z=\pi/L_z$)

Нормировка на $P_{src}=-\frac12\mathrm{Re}\int j^{ext*}_yE_y\,dx$. Квадратура Гаусса–Лежандра до $40/\mathrm{Im}K$.

| $\omega$, рад/с | объём правой стенки $\int_a^\infty f^M_x$ | скачок $[T_{xx}]_a$ | $F_r=-T^{vac}_{xx}(a)$ | $F_l=T^{vac}_{xx}(-a)$ | $F_l+F_r$ | $P_{src}$, эрг/(с·см²) | $\frac{P_{src}-P_{diss}}{P_{src}}$ | невязка «объём = $-T^M(a^+)$» | невязка Э–Л |
|---|---|---|---|---|---|---|---|---|---|
| 1.0000e+06 | 4.4712e-06 | -4.6881e-06 | -2.1684e-07 | 2.1684e-07 | -1.0e-17 | 5.3151e-09 | 2e-18 | 1e-16 | 1e-16 |
| 2.0000e+06 | 6.1553e-06 | -6.1858e-06 | -3.0446e-08 | 3.0446e-08 | -7.3e-18 | 3.7582e-09 | 2e-18 | 1e-16 | 1e-16 |
| 3.0000e+06 | 7.2270e-06 | -7.2044e-06 | 2.2585e-08 | -2.2585e-08 | -5.9e-18 | 3.0685e-09 | 2e-18 | 1e-16 | 1e-16 |
| 4.0000e+06 | 7.9035e-06 | -7.8565e-06 | 4.7083e-08 | -4.7083e-08 | -5.1e-18 | 2.6573e-09 | 2e-18 | 1e-16 | 1e-16 |
| 5.0000e+06 | 8.2960e-06 | -8.2356e-06 | 6.0352e-08 | -6.0352e-08 | -4.6e-18 | 2.3768e-09 | 2e-18 | 1e-16 | 1e-16 |
| 6.0000e+06 | 8.4835e-06 | -8.4157e-06 | 6.7841e-08 | -6.7841e-08 | -4.2e-18 | 2.1696e-09 | 2e-18 | 1e-16 | 1e-16 |
| 7.0000e+06 | 8.5269e-06 | -8.4550e-06 | 7.1947e-08 | -7.1947e-08 | -3.9e-18 | 2.0087e-09 | 2e-18 | 1e-16 | 1e-16 |
| 8.0000e+06 | 8.4719e-06 | -8.3980e-06 | 7.3946e-08 | -7.3946e-08 | -3.6e-18 | 1.8789e-09 | 2e-18 | 1e-16 | 1e-16 |
| 9.0000e+06 | 8.3523e-06 | -8.2777e-06 | 7.4600e-08 | -7.4600e-08 | -3.4e-18 | 1.7715e-09 | 2e-18 | 1e-16 | 1e-16 |
| 1.0000e+07 | 8.1919e-06 | -8.1175e-06 | 7.4387e-08 | -7.4387e-08 | -3.2e-18 | 1.6805e-09 | 2e-18 | 1e-16 | 1e-16 |
| 1.7988e+11 | 1.8650e-07 | -1.8441e-07 | 2.0863e-09 | -2.0863e-09 | -2.4e-20 | 1.2529e-11 | 2e-18 | 1e-16 | 1e-16 |

**Вывод из таблицы.** Объём правой стенки ($\sim5\cdot10^{-6}$ Н/кВт) почти полностью компенсирует скачок на поверхности. Полная сила на стенку $F_r=-T^{vac}_{xx}(a)\sim10^{-8}$–$10^{-7}$ Н/кВт.
Левая стенка вместе с источником получает $-F_r$, сумма — ноль. Сравнение: фотонная ракета даёт $1/c=3.3\cdot10^{-6}$ Н/кВт.

---

## 4. Найденные ошибки (номера ячеек `MenDrive_qnm.ipynb`)

| # | где | ошибка | влияние |
|---|---|---|---|
| E1 | `mendrivepylib/variables/common.py` (нет `assume(mu_r_*_c,'complex')`), ячейки 291, 293 (`.simplify()`) | $\overline{\mu'_r}\to\mu'_r$ в $F_{P\times\mathrm{rot}E,r}$ и $F_{I\times\mathrm{rot}B,r}$ | **+3466 Н/кВт и $1/\omega^2$** (знак и величина $\sim\mu''^2$) |
| E2 | ячейка 396 (`force_dyn_left/right`) | сумма 12 членов разных формулировок: Абрагам + Эйнштейн–Лауб + Смайт + скачок Минковского | двойной и тройной счёт; без E1 было бы $-3517$ Н/кВт |
| E3 | ячейки 291–319 | $\mathbf P'=\frac{\varepsilon'-1}{4\pi}\mathbf E$, $\mathbf M'=\frac{\mu'-1}{4\pi}\mathbf H$ с $\sigma$ внутри | ток проводимости и магнитный ток посчитаны как поляризация и намагниченность ($\propto\sigma/\omega$) |
| E4 | ячейка 396 (`svd_X4`), 476 (`scan_omega_re`) | при вещественной $\omega$ берётся ядро однородной $M_4$; $B_{ext}j_0$ не используется | поля не удовлетворяют сшивке $E_y$; нормировка $\lVert X\rVert=1$ произвольна |
| E5 | ячейки 21, 47, 75, 81–115, 232, 473, 477 | источник $j_0e^{-iK_l(x+a)}$ резонансен, а $H_z$ оставлен однородным: нет секулярного члена, в `eq_Ey_l` стоит $4\pi$ вместо $2\pi$ | амплитуды вынужденного решения $\times2$; в левой стенке нарушен закон Фарадея |
| E6 | ячейка 278 (`E_y_l`, `H_x_l` через `M_left`) | член с `je0_ext` есть и в столбце `'Ez'`, и в столбце `'Hz'`, поэтому удвоен ($8\pi$) | силы левой стенки (относительно $10^{-5}$) |
| E7 | ячейки 392, 396 | $W$ — потери только в слоях $h_l,h_r$, а силы интегрируются до $\infty$; $P_{src}$ не считается | при $\omega=10^6$ учтено 12 % потерь ($\delta_l=0.016$ см, $\delta_r=0.094$ см $\gg h$) |
| E8 (латентная) | ячейки 21, 280 (`jm_ext_z_profile`) | правый профиль $e^{-iK_r(x+a)}$ вместо $e^{-iK_r(a-x)}$: $\int_a^\infty$ расходится | в `driving_mode_1` не активна ($F_{jm,ext,D}=0$) |
| E9 (замечание) | ячейки 551–760 (QNM) | при комплексной $\omega=\omega'+i\omega''$ $\langle\partial_t\mathbf g\rangle=2\omega''\langle\mathbf g\rangle\ne0$ | баланс «объём = тензор» на QNM-частоте содержит член $2\omega''\int\langle g_x\rangle dx$; это не стационарная тяга |

Проверено, **не является ошибкой**:
* `sign_K`: `fix_K_sign` выбирает $\mathrm{Im}K>0$ при `sign_K_1_subs`, а анзацы $(x+a)$ слева и $(x-a)$ справа затухают.
* $\mathbf j_e\times\mathbf H$ (`p_Ht_jeH`) в сумму не входит.
* Перекрёстные члены $\sigma_e\sigma_m/\omega$ в $F_{jmD}+p_{Bt,jeB}$ сокращаются (§2).

---

## 5. Исправления (код Sage для вставки)

**F1 (E1).** Новая ячейка сразу после ячейки 10 (`from mendrivepylib.variables.qnm import *`). После неё перевыполнить ячейки 268–347:
```python
# FIX E1: без этого Maxima в .simplify() считает mu_r_*_c вещественными и выбрасывает conjugate()
for _v in [mu_r_xx_c, mu_r_yy_c, mu_r_zz_c, mu_r_yz_c, mu_r_zy_c, mu_l_yz_c, mu_l_zy_c,
           epsilon_r_xy_c, epsilon_r_yx_c, epsilon_r_xz_c, epsilon_r_zx_c, epsilon_r_yz_c, epsilon_r_zy_c]:
    assume(_v, 'complex')
# контроль (после ячейки 293): в правой части не должно быть mu_r_zz_c^2
assert 'conjugate(mu_r_zz_c)' in str(f_I_rotB_conj_r) and 'conjugate(mu_r_xx_c)' in str(f_I_rotB_conj_r)
```
Альтернатива без правки библиотеки: в ячейках 291 и 293 убрать `.simplify()`.

**F2 (E2, E3).** В `calc_thrust` (ячейка 396) заменить определения `force_dyn_left` / `force_dyn_right`:
```python
# FIX E2/E3: Минковский -- объём только на свободные токи + скачок тензора на границе с вакуумом
force_dyn_right = (F_jmD_conj_right + p_Bt_jeB_conj_right + F_jm_ext_D_conj_right
                   + f_jump_maxwell_right).n(prec=128)
force_dyn_left  = (F_jmD_conj_left + p_Bt_jeB_conj_left + p_Bt_je_ext_B_conj_left
                   + f_jump_maxwell_left).n(prec=128)
# контроль: объём = тензор (тождество Минковского) и полная сила = 0
_tol = 1e-8*abs(f_maxwell_right_wall) + 1e-40
if abs(force_dyn_right - f_maxwell_right_wall) > _tol or abs(force_dyn_left - f_maxwell_left_wall) > 1e-6*abs(f_maxwell_left_wall) + 1e-40:
    print('FAIL: объём != тензор', force_dyn_right, f_maxwell_right_wall, force_dyn_left, f_maxwell_left_wall)
# Абрагамовы/динамические/смайтовские члены оставить только как диагностику (в result), в сумму НЕ включать
```

**F3 (E4).** В `calc_thrust` вместо `X4, ans = svd_X4(M4_d_cur)`:
```python
_je0 = [e.rhs() for e in omega_sol_cur if e.lhs() == je0_ext]
_je0 = complex128(_je0[0]) if _je0 else 0
if _je0 != 0:      # вынужденный режим (вещественная omega): M4 X = B_ext*je0
    B_d = vector([b.subs(K_v_d).subs(K_l_d).subs(K_r_d).subs(kappa_vacuum_sol).subs(epsilon_mu_subs)
                    .subs(k_subs).subs(dv).subs(omega_sol_cur).subs(sign_K_1_subs) for b in B4_K])
    M_n = matrix(ComplexField(128), M4_d_cur.n(prec=128)); B_n = vector(ComplexField(128), B_d.n(prec=128))
    X4 = M_n.solve_right(B_n)                    # после F4 B4_K[2] содержит 2*pi
    ans = {'residual': float((M_n*X4 - B_n).norm()/B_n.norm())}
else:              # QNM: ядро однородной системы
    X4, ans = svd_X4(M4_d_cur)
```

**F4 (E5).** Ячейка 47, ветка TE с `driving_mode_1`. Секулярный член (получен из $H_z''+K_l^2H_z=-\frac{4\pi}{c}\partial_xj_y$):
```python
if is_TE:
    if driving_mode_1 == driving_mode and not get_use_phase_y():
        H_left_conductor_sol_z (x) = (B_H_zl - 2*pi*je0_ext/c*(x + a)) \
                                     * exp(-I * sign_K_H_l * K_H_left_conductor*(x+a))
    else:
        H_left_conductor_sol_z (x) = (B_H_zl*exp(-I * sign_K_H_l * K_H_left_conductor*(x+a)))
```
После этого ячейки 102–131 дадут
$E_y=\frac{c\mu'\omega K_l(B_{H,zl}-\frac{2\pi j_0}{c}(x+a))+2\pi i\mu'\omega j_0}{c^2k_z^2-\varepsilon'\mu'\omega^2}e^{-iK_l(x+a)}$,
а `GenerateRHS` — $B_{ext}[2]=-\frac{2\pi i\mu'\omega}{D}$. Проверка: $X(\text{F4})=X(\text{ячейка 477})/2$ до $10^{-33}$.
Шаблоны `.collect(e^(...))` и `m_integrate` (ячейка 298) рассчитаны на чистые экспоненты. Для членов с $(x+a)$ интегрировать напрямую:
`integrate(f.expand(), x, -Infinity, -a)` или численно `numerical_integral(real(f_num), -A-40/ImK, -A)`.
Минимальная замена без перевывода: в ячейке 476 после `B_ext_final = ...` поставить `B_ext_final = B_ext_final/2`. Это делает точными поля в зазоре и в правой стенке; поля внутри левой стенки остаются без секулярного члена.

**F5 (E6).** Ячейка 278: для TE ($B_{E,zl}=0$) столбец `'Ez'` не добавлять:
```python
E_y_l(x) = (M_left[('Ey','Ez')]*B_E_zl*exp(-I*K_E_l*(a+x))*phase_y_TM if is_TM else 0) \
         +  M_left[('Ey','Hz')]*B_H_zl*exp(-I*K_H_l*(a+x))*phase_y_TE
H_x_l(x) = (M_left[('Hx','Ez')]*B_E_zl*exp(-I*K_E_l*(a+x))*phase_y_TM if is_TM else 0) \
         +  M_left[('Hx','Hz')]*B_H_zl*exp(-I*K_H_l*(a+x))*phase_y_TE
```

**F6 (E7).** В `calc_thrust` (ячейка 396) считать потери по всей толщине скин-слоя. В линейном стационарном режиме это равно $P_{src}$:
```python
# FIX E7: вместо h_conductor_l/r (слои 1.9e-4 и 6.3e-3 см) -- полубесконечные стенки
ImKl = float(K_l_d[0].rhs().imag()); ImKr = float(K_r_d[0].rhs().imag())
h_l = 40/ImKl          # exp(-80) -- полное затухание
h_r = 40/ImKr
# (далее как было: W_r_val = Wr_int(A+h_r) - Wr_int(A), W_l_val = -Wl_int(-A-h_l) + Wl_int(-A))
# контроль мощности источника (после F4): P_src = -1/2 Re int je_ext^* E_y dx  ==  W_total_ergs
```
Н/кВт $=10^5\,F/P_{src}$ (как в ноутбуке: `force_N/W_total_W*1e3`).

**F7 (E8).** Ячейки 21 и 280: `jm_ext_z_profile(x) = jm0_ext * exp(-I*sign_K_H_r*K_H_right_conductor*(a - x)) * phase_y_TE`.

---

## 6. Задача B: согласованный конфиг C++ (`mendrive_llg/MenDriveCpp/standalone`)

### 6.1 Единицы и модель C++ (по `mendrive_core.cpp`, `main.cpp`)

* Сетка Йи на $x\in[-(a+h_l),\,a+h_r]$, $\Delta x=(2a+h_l+h_r)/N$, $\Delta t=$ `dt_frac`$\cdot\Delta x$; на внешних концах $E=0$ (PEC). Вакуум $|x|\le a$.
* $c=1$, время $\tau=ct$. Левая стенка: `eps_left`, `sigma_e_left`, $\mu=1$, $\sigma_m=0$. Феррит (Linear): `eps_ferrite`, `sigma_e_ferrite`, `mu_zz_lin` (для $H_z$), `mu_yy_lin`, $\sigma_m=$ `sigma_m_leak`.
  Полунеявные коэффициенты $(\varepsilon\mp2\pi\sigma\Delta t)$ и $(\mu\mp2\pi\sigma_m\Delta t)$ означают $\varepsilon'=\varepsilon+4\pi i\sigma/\omega$, $\mu'=\mu+4\pi i\sigma_m/\omega$.
* Источник `je_y_left`: $j=\text{amp}\,e^{-\xi/\delta}\sin(\omega\tau-k\xi)$, $\xi=-a-x$, то есть $\tilde j=i\,\text{amp}\,e^{-iK_s(x+a)}$, $K_s=k+i/\delta$.
* Тяга `dTxx` $=\langle p(a)\rangle-\langle p(-a)\rangle$, $p=(E^2+H^2)/8\pi$ (`Txx_at`, давление вакуума). Сила на одну стенку — столбец `Txx` в `fields_harm_n0.csv`.

**Отображение (масштаб $L_0=1$ см, $T_0=L_0/c$).** $x_{code}=x/L_0$, $\omega_{code}=\omega L_0/c$, $\sigma_{code}=\sigma_{СГС}L_0/c$ (и $\sigma_e$, и $\sigma_m$), $\varepsilon,\mu$ без изменений, поля гауссовы.
$j_{code}=j\,L_0/c$, $P_{code}=P_{СГС}/c$, откуда
$$\text{Н/кВт}=10^5\frac{F}{P}=10^5\frac{\texttt{dTxx}}{c\,P_{src}},$$
что совпадает с `force_per_kW_src` при `length_unit_cm=1`, `time_unit_s=1/c`. **По умолчанию** (`length_unit_cm=1`, `time_unit_s=1`) `force_per_kW*` завышены в $c=3\cdot10^{10}$ раз.
Значения из `log.txt` (+20.3 при $N=80$ и $-1.84$ Н/кВт при $N=160$; Preisach, magnetic_right) физически равны $+6.8\cdot10^{-10}$ и $-6.1\cdot10^{-11}$ Н/кВт.
Смена знака при измельчении сетки показывает ошибку дискретизации: точное значение равно 0 (§2).

| ноутбук (СГС) | C++ | значение |
|---|---|---|
| $a=0.5$, $h_l$, $h_r$ | `a`, `h_l`, `h_r` | 0.5, 0.0001875, 0.00625 |
| $\varepsilon_l=1$, $\sigma_{e,l}=5.625\cdot10^{17}$ | `eps_left`, `sigma_e_left` | 1, 18762980.354896054 |
| $\varepsilon_r=16$, $\sigma_{e,r}=9.9\cdot10^{12}$ | `eps_ferrite`, `sigma_e_ferrite` | 16, 330.22845424617054 |
| $\mu_r=100$, $\sigma_{m,r}=6.48\cdot10^{7}$ | `ferrite_model=Linear`, `mu_zz_lin`, `mu_yy_lin`, `sigma_m_leak` | 100, 100, 0.0021614953368840253 |
| $\mathrm{Re}\,\omega_{QNM}=1.7988142\cdot10^{11}$ | `omega0` | 6.000198342077812 |
| $j^{ext}_y=j_0e^{-iK_l(x+a)}$ | `excitation_mode=je_y_left`, `skin_depth=1/Im K_l`, `src_k_phase=Re K_l`, `src_sigma_fix=1` | 3.75990053432236e-05, 26596.448930483315 |
| $k_z=\pi/L_z$ | — (в 1D $k_z=0$) | вклад $(ck_z/\omega)^2\approx3\cdot10^{-4}$ |

`src_sigma_fix=1` обязателен: при $2\pi\sigma_l\Delta t\approx450$ без деления на $(\varepsilon+2\pi\sigma\Delta t)$ источник был бы завышен примерно в 450 раз.

### 6.2 Что C++ не умеет и минимальные расширения (только описание)

1. **$k_z\ne0$** (1D-модель). Расширение: член $-k_z^2$ в уравнение для $H_x$ ($H_x=-ck_zE_y/(\omega\mu)$) — в $T_{xx}$ добавится $|H_x|^2$. Для $\omega\approx\omega_{QNM}$ не нужно.
2. **Полубесконечные стенки.** Стенки конечные, за ними PEC. При $h_l=5.0\delta_l$, $h_r=6.9\delta_r$ разница в зазоре $7\cdot10^{-4}$ (проверено). Для $\omega\le10^7$ это неприменимо ($\delta\gg h$).
3. **Равномерная сетка.** $\delta_l=3.76\cdot10^{-5}$ см требует $N\gtrsim6.5\cdot10^4$; основной конфиг — $N=131072$ ($\delta_l/\Delta x=4.9$).
   Расширение: неравномерная сетка (сгущение у $x=-a$) или граница Леонтовича $E_y=-Z H_z$, $Z=(1-i)\sqrt{\omega/8\pi\sigma_l}$ для левого металла.
4. **Сила на отдельную стенку** в `summary.txt` не выводится. Расширение: печатать $\langle p(a)\rangle$, $\langle p(-a)\rangle$ отдельно. Сейчас её даёт `fields_harm_n0.csv` (n=0 не нормируется).
5. `Txx_at` смешивает $E^n$ и $H^{n+1/2}$ (ошибка $O(\Delta t)$) и усредняет $H$ по соседним ячейкам (ошибка $O(\Delta x)$). Расширение: центрировать $H$ по времени.
6. Частоты скана ячейки 479 ($\omega=10^6..10^7$, $\omega_{code}=3.3\cdot10^{-5}..3.3\cdot10^{-4}$) во FDTD недостижимы: $\sim10^{10}$–$10^{11}$ шагов на период.
7. `linear1d.solve_linear` **переполняется** при физических скин-слоях: `OverflowError: math range error` в строке `C = j0*cmath.exp(-1j*Ks*a)`, где $|e^{-iK_sa}|=e^{13298}$.
   Минимальная правка: `C = j0` и `part = up*exp(q*(x - xref))`, $x_{ref}=-a$ для левой стенки и $+a$ для правой.
   С этой правкой `linear1d` (halfspace, $K_s=K_l(1+10^{-7})$) совпадает с точным решением ноутбука: при $k_z=0$ до $5\cdot10^{-8}$, с $k_z=\pi/L_z$ — до $3\cdot10^{-4}$ (по $E_y$; по $H_z$ до $10^{-6}$).
   Пока правки нет, вместо `python3 linear1d.py out_dir` используйте `export_qnm_profiles_sage.py --kz0`.

### 6.3 Созданные файлы

* `examples/qnm_drive_je_y_left_linear.cfg` — основной, $N=131072$, $\omega_{code}=6.000198$, 60 периодов (запись с 40-го), `fields_harm=0,1`, `fields_harm_normalize=source`, `fields_snap_phases=0,0.25,0.5,0.75`.
* `examples/qnm_drive_je_y_left_linear_N80.cfg`, `_N160.cfg`, `_N320.cfg` — те же параметры. Это **контрольные** прогоны: $\Delta x>h_l$, в левой стенке нет внутренних узлов, источник не действует ($P_{src}=0$, `nan`).
* `qnm_sweep.sh` — сетка `N_LIST` × `W_LIST`. Для каждой $\omega$ пересчитывает `skin_depth` и `src_k_phase`, умеет `JOBS`, `DRY=1`, `SUMMARY_ONLY=1`, пропускает готовые прогоны.
  Пишет `out/qnm_sweep/summary_table.{csv,md}`: `dTxx, P_src, P_diss, P_num, rel_src, NkW_thrust, NkW_right, NkW_right_exact, P_src_ratio`.
* `export_qnm_profiles_sage.py` — профили $\tilde E_y,\tilde H_z$ в CSV для `plot.py --compare`: (1) из ноутбука с вынужденным решением, (2) из готовых `Ey_*_d/Hz_*_d`, (3) без Sage — точное решение. Проверено: совпадает с независимым решателем до $10^{-15}$.

Ожидаемые значения для C++ (точное 1D, $k_z=0$, amp=1):

| конфигурация | $\omega_{code}$ | `skin_depth` | `src_k_phase` | $P_{src}$ (код, amp=1) | $\langle p\rangle(a)$ (код) | Н/кВт на правую стенку | Н/кВт суммарно |
|---|---|---|---|---|---|---|---|
| QNM0_Re ($\omega=8.82769e+10$) | 2.944600605 | 5.367177e-05 | 18631.767122 | 5.36196e-13 | 4.36234e-16 | 2.7138e-09 | 0 |
| 0.9*QNM1 ($\omega=1.61893e+11$) | 5.400178508 | 3.963283e-05 | 25231.606856 | 3.96069e-13 | 3.97359e-17 | 3.3465e-10 | 0 |
| QNM1_Re ($\omega=1.79881e+11$) | 6.000198342 | 3.759901e-05 | 26596.448930 | 3.75620e-13 | 2.34937e-16 | 2.0863e-09 | 0 |
| 1.1*QNM1 ($\omega=1.97870e+11$) | 6.600218176 | 3.584925e-05 | 27894.591004 | 3.58192e-13 | 4.03485e-17 | 3.7574e-10 | 0 |
| QNM2_Re ($\omega=2.72016e+11$) | 9.073479329 | 3.057541e-05 | 32706.026292 | 3.05450e-13 | 1.69899e-16 | 1.8554e-09 | 0 |

Критерии сходимости по $N$: `P_src_ratio`$\to1$, `rel_dif_src`$\to0$, `NkW_right`$\to$`NkW_right_exact`, `NkW_thrust`$\to0$ (как $O(\Delta x)$).

### 6.4 Команды в docker

```sh
cd /mnt/raid/articles/sagemath_docker_build/work/mendrive_llg/MenDriveCpp/standalone
python3 -m pip install --user numpy matplotlib      # нужно для make test / plot.py (в log.txt: No module named 'numpy')
#   или без pip: sage -python tests/compare_ctypes.py ; sage -python plot.py ...  (Makefile не меняем, PYTHON там нет)
make CXX=g++
make test                                            # если python3 без numpy: sage -python tests/compare_ctypes.py && sage -python tests/test_fields_linear.py
./mendrive -c examples/qnm_drive_je_y_left_linear_N160.cfg   # секунды; контроль: P_src=0 -> nan
./mendrive -c examples/qnm_drive_je_y_left_linear.cfg        # основной, ~1-2.5 ч
python3 export_qnm_profiles_sage.py --kz0 --out out/qnm_drive_je_y_left_linear/nb_exact_kz0.csv
python3 export_qnm_profiles_sage.py --kz0 --notebook-source --out out/qnm_drive_je_y_left_linear/nb_as_notebook_kz0.csv
python3 plot.py out/qnm_drive_je_y_left_linear --compare out/qnm_drive_je_y_left_linear/nb_exact_kz0.csv --fit none --norm none
N_LIST="16384 65536 131072" W_LIST="6.000198342077812 9.073479328679792" JOBS=2 bash qnm_sweep.sh
bash qnm_sweep.sh                                    # полная сетка по умолчанию (долго); DRY=1 -- только план
cd /mnt/raid/articles/sagemath_docker_build/work/rezonatory_i_volnovody
jupyter nbconvert --to notebook --execute check_force_tensor_vs_volume.ipynb --output check_force_tensor_vs_volume.out.ipynb
#   (или sage -n jupyter и Run All; первая ячейка выключает препарсер Sage)
```
Из ноутбука `MenDrive_qnm.ipynb` (после F1–F5), в конце:
```python
import sys; sys.path.insert(0, '/mnt/raid/articles/sagemath_docker_build/work/mendrive_llg/MenDriveCpp/standalone')
import export_qnm_profiles_sage as E
F = E.driven_fields_from_notebook(globals(), omega_drive=1.7988142094590323e11, je0=1e6, source_fix=True)
E.export_notebook_fields(globals(), F, '/mnt/raid/articles/sagemath_docker_build/work/mendrive_llg/MenDriveCpp/standalone/out/qnm_drive_je_y_left_linear/nb_profile.csv', je0=1e6)
```
(`source_fix=True` нужен, пока не сделана F4.)

---

## ВЫВОД

* Расхождение «объём 3466 Н/кВт против тензора $-0.0071$» вызвано не физикой, а ошибками в ноутбуке.
  Главная — потерянное $\overline{\mu'_r}$ в $F_{I\times\mathrm{rot}B}$ (E1): член $\propto(4\pi\sigma_m/\omega)^2|H_z|^2/16\pi$, отсюда $1/\omega^2$.
  Вторая — сумма несовместимых формулировок с двойным счётом токов проводимости (E2, E3).
* При правильной силе объём + поверхность = тензор, тождество выполнено до $10^{-16}$. Полная сила на систему равна нулю. Сила на каждую стенку $\lesssim10^{-7}$ Н/кВт, левая и правая противоположны.
* Вынужденные поля ноутбука надо пересчитать: решать $M_4X=B_{ext}j_0$ вместо SVD, с источником $2\pi$ и секулярным членом. Нормировать на $P_{src}$.

## ОТКРЫТО

1. После F1–F6 перезапустить ячейку 479 и убедиться, что `thrust_N_per_kW` равна 0 до округления, а `thrust_N_per_kW_right` = `N_per_kW_maxwell_right_wall`.
2. QNM-ячейки 551–760: при комплексной $\omega$ в балансе нужен член $\langle\partial_t g\rangle=2\omega''\langle g\rangle$. Если смотреть «тягу» QNM, нормировать её на энергию моды, а не на мощность.
3. Прогон C++ на $N=65536$ и $131072$. Ожидание: `NkW_thrust`$\to0$, `NkW_right`$\to2.09\cdot10^{-9}$ при $\omega_{code}=6.0002$. Отклонение покажет порядок дискретизационной «тяги» схемы (как в `log.txt`).
4. Для $\omega\le10^7$ (скан ячейки 479) C++ неприменим. Сравнивать только с линейной моделью (`check_force_tensor_vs_volume.ipynb`, §6).
5. Нелинейные режимы (Preisach/JA): $\langle\texttt{dTxx}\rangle=0$ для любого периодического режима (вакуум между стенками). Ненулевое значение — признак нестационарности или дискретизации ($P_{num}$, `Txx_at`).
