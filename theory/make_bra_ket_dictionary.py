# -*- coding: utf-8 -*-
"""
Генератор PDF: словарь бра-кет <-> интегралы Ландау-Давыдова
+ вывод коэффициента Боголюбова из гауссова интеграла
+ формулы Степановского (4)-(11) на языке волновых функций.

Без LaTeX: только ReportLab (pip install reportlab).
Нужен TTF-шрифт с кириллицей и символами (DejaVuSans): скрипт ищет его сам,
либо укажите путь в переменной FONT_DIR.
Совместимо с Python 3.5+ (без f-строк).
Запуск:  python3 make_bra_ket_dictionary.py  [выходной.pdf]
"""
import os, sys, glob
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, white

FONT_DIR = None  # можно задать вручную, например "/Library/Fonts/"

def find_font_dir():
    if FONT_DIR:
        return FONT_DIR
    pats = ["/usr/share/fonts/**/DejaVuSans.ttf",
            "/Library/Fonts/**/DejaVuSans.ttf",
            os.path.expanduser("~/Library/Fonts/**/DejaVuSans.ttf")]
    try:
        import matplotlib
        pats.insert(0, os.path.join(os.path.dirname(matplotlib.__file__),
                                    "mpl-data", "fonts", "ttf", "DejaVuSans.ttf"))
    except ImportError:
        pass
    for p in pats:
        hits = glob.glob(p, recursive=True)
        if hits:
            return os.path.dirname(hits[0]) + os.sep
    sys.exit("Не найден DejaVuSans.ttf: укажите FONT_DIR в начале скрипта")

D = find_font_dir()
for name, f in [("DV", "DejaVuSans.ttf"), ("DVB", "DejaVuSans-Bold.ttf"),
                ("DVS", "DejaVuSerif-Bold.ttf")]:
    pdfmetrics.registerFont(TTFont(name, D + f))

P = HexColor("#1E2761"); ACC = HexColor("#8A3B12")
def st(**k): return ParagraphStyle("x", **k)
T   = st(fontName="DVS", fontSize=17, leading=22, textColor=P, spaceAfter=4)
Sub = st(fontName="DV", fontSize=9.5, leading=13, textColor=HexColor("#555555"), spaceAfter=10)
H   = st(fontName="DVS", fontSize=12, leading=16, textColor=P, spaceBefore=10, spaceAfter=5)
C   = st(fontName="DV", fontSize=8.6, leading=11.5)
CB  = st(fontName="DVB", fontSize=8.6, leading=11.5, textColor=white)
N   = st(fontName="DV", fontSize=9.3, leading=13.5, spaceAfter=4)
EQ  = st(fontName="DV", fontSize=9.8, leading=14, leftIndent=14, spaceBefore=2, spaceAfter=5,
         textColor=HexColor("#111111"))
BOX = st(fontName="DV", fontSize=9.3, leading=13.5, textColor=ACC, spaceBefore=4, spaceAfter=6)

def tbl(rows, head=("Бра-кет", "Язык Ландау–Давыдова", "Смысл")):
    d = [[Paragraph(x, CB) for x in head]]
    d += [[Paragraph(a, C), Paragraph(b, C), Paragraph(c, C)] for a, b, c in rows]
    t = Table(d, colWidths=[52*mm, 72*mm, 52*mm], repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), P),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [white, HexColor("#F3F4F8")]),
        ("GRID", (0, 0), (-1, -1), 0.4, HexColor("#C8CBD6")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4)]))
    return t

SEC = [
("1. Состояния и скалярное произведение", [
("∣ψ⟩", "ψ(x)", "состояние (волновая функция)"),
("⟨ψ∣", "ψ*(x)", "сопряжённое состояние"),
("⟨x∣ψ⟩", "ψ(x)", "волновая функция в x-представлении"),
("⟨p∣ψ⟩", "φ(p) = (2πħ)<sup>−1/2</sup>∫e<sup>−ipx/ħ</sup>ψ(x)dx", "импульсное представление"),
("⟨φ∣ψ⟩", "∫φ*(x)ψ(x)dx", "скалярное произведение (интеграл перекрытия)"),
("⟨ψ∣ψ⟩ = 1", "∫∣ψ∣²dx = 1", "нормировка"),
("⟨φ∣ψ⟩* = ⟨ψ∣φ⟩", "(∫φ*ψdx)* = ∫ψ*φdx", "свойство сопряжения")]),
("2. Базисы, ортогональность, полнота", [
("⟨m∣n⟩ = δ<sub>mn</sub>", "∫ψ<sub>m</sub>*ψ<sub>n</sub>dx = δ<sub>mn</sub>", "ортонормированность"),
("⟨x∣x′⟩ = δ(x−x′)", "—", "собственные функции координаты"),
("⟨p∣p′⟩ = δ(p−p′)", "∫ψ<sub>p</sub>*ψ<sub>p′</sub>dx = δ(p−p′)", "нормировка непрерывного спектра"),
("⟨x∣p⟩", "(2πħ)<sup>−1/2</sup>e<sup>ipx/ħ</sup> = ψ<sub>p</sub>(x)", "плоская волна"),
("Σ<sub>n</sub>∣n⟩⟨n∣ = 1", "Σ<sub>n</sub>ψ<sub>n</sub>(x)ψ<sub>n</sub>*(x′) = δ(x−x′)", "условие полноты"),
("∫∣x⟩⟨x∣dx = 1", "—", "полнота координатного базиса"),
("c<sub>n</sub> = ⟨n∣ψ⟩", "c<sub>n</sub> = ∫ψ<sub>n</sub>*ψdx", "коэффициент разложения"),
("∣ψ⟩ = Σc<sub>n</sub>∣n⟩", "ψ(x) = Σc<sub>n</sub>ψ<sub>n</sub>(x)", "разложение по собственным функциям")]),
("3. Операторы", [
("Â∣ψ⟩", "Âψ(x)", "действие оператора"),
("⟨x∣x̂∣ψ⟩", "xψ(x)", "оператор координаты"),
("⟨x∣p̂∣ψ⟩", "−iħ dψ/dx", "оператор импульса"),
("⟨x∣Â∣x′⟩", "A(x,x′):  Âψ(x) = ∫A(x,x′)ψ(x′)dx′", "интегральное ядро"),
("⟨m∣Â∣n⟩ = A<sub>mn</sub>", "∫ψ<sub>m</sub>*Âψ<sub>n</sub>dx", "матричный элемент"),
("⟨ψ∣Â∣ψ⟩ = ⟨A⟩", "∫ψ*(x)Âψ(x)dx", "среднее значение оператора"),
("Â∣n⟩ = a<sub>n</sub>∣n⟩", "Âψ<sub>n</sub> = a<sub>n</sub>ψ<sub>n</sub>", "задача на собственные значения"),
("⟨φ∣Â<sup>†</sup>∣ψ⟩ = ⟨ψ∣Â∣φ⟩*", "∫φ*Â<sup>†</sup>ψdx = (∫ψ*Âφdx)*", "эрмитово сопряжение"),
("Â = Â<sup>†</sup>", "∫φ*Âψdx = ∫(Âφ)*ψdx", "эрмитовость"),
("[Â,B̂] = ÂB̂ − B̂Â", "то же на функциях; [x̂,p̂] = iħ", "коммутатор"),
("Tr Â = Σ⟨n∣Â∣n⟩", "ΣA<sub>nn</sub> = ∫A(x,x)dx", "след")]),
("4. Кет-бра: проекторы и матрица плотности", [
("∣ψ⟩⟨φ∣", "ядро ψ(x)φ*(x′)", "оператор, а не число"),
("∣n⟩⟨n∣", "ψ<sub>n</sub>(x)ψ<sub>n</sub>*(x′)", "проектор на состояние n"),
("ρ̂ = ∣ψ⟩⟨ψ∣", "ρ(x,x′) = ψ(x)ψ*(x′)", "матрица плотности чистого состояния (ЛЛ §14)"),
("ρ̂ = Σw<sub>k</sub>∣ψ<sub>k</sub>⟩⟨ψ<sub>k</sub>∣", "Σw<sub>k</sub>ψ<sub>k</sub>(x)ψ<sub>k</sub>*(x′)", "смешанное состояние"),
("⟨A⟩ = Tr(ρ̂Â)", "∫[Âρ(x,x′)]<sub>x′=x</sub>dx", "среднее через матрицу плотности")]),
("5. Время и многочастичные системы", [
("iħ∂<sub>t</sub>∣ψ⟩ = Ĥ∣ψ⟩", "iħ∂ψ(x,t)/∂t = Ĥψ(x,t)", "уравнение Шрёдингера"),
("∣ψ(t)⟩ = e<sup>−iĤt/ħ</sup>∣ψ(0)⟩", "ψ(x,t) = Σc<sub>n</sub>ψ<sub>n</sub>(x)e<sup>−iE<sub>n</sub>t/ħ</sup>", "эволюция"),
("Â(t) = e<sup>iĤt/ħ</sup>Âe<sup>−iĤt/ħ</sup>", "—", "картина Гейзенберга (ЛЛ §13)"),
("∣a,b⟩ = ∣a⟩⊗∣b⟩", "ψ<sub>a</sub>(x₁)ψ<sub>b</sub>(x₂)", "система двух частиц"),
("⟨a,b∣c,d⟩", "∬ψ<sub>a</sub>*(x₁)ψ<sub>b</sub>*(x₂)ψ<sub>c</sub>(x₁)ψ<sub>d</sub>(x₂)dx₁dx₂", "двойной интеграл"),
("Tr<sub>2</sub>ρ̂", "ρ(x₁,x₁′) = ∫ρ(x₁,x₂;x₁′,x₂)dx₂", "частичный след")]),
("6. Осциллятор, операторы a, a† и Боголюбов", [
("∣n⟩", "ψ<sub>n</sub>(x) = N<sub>n</sub>H<sub>n</sub>(ξ)e<sup>−ξ²/2</sup>,  ξ = x√(mω/ħ)", "n-й уровень = n квантов"),
("∣0⟩", "ψ<sub>0</sub>(x) = (mω/πħ)<sup>1/4</sup>e<sup>−mωx²/2ħ</sup>", "вакуум = гауссиана"),
("â", "√(mω/2ħ)·(x + (ħ/mω)d/dx)", "понижающий оператор"),
("â<sup>†</sup>", "√(mω/2ħ)·(x − (ħ/mω)d/dx)", "повышающий оператор"),
("â∣n⟩ = √n∣n−1⟩", "âψ<sub>n</sub> = √n ψ<sub>n−1</sub>", "понижение уровня"),
("â∣0⟩ = 0", "âψ<sub>0</sub> = 0 (дифференцированием гауссианы)", "ниже вакуума уровня нет"),
("⟨n∣â<sup>†</sup>â∣n⟩ = n", "∫ψ<sub>n</sub>*(â<sup>†</sup>âψ<sub>n</sub>)dx = n", "оператор числа квантов"),
("[â,â<sup>†</sup>] = 1", "â(â<sup>†</sup>ψ) − â<sup>†</sup>(âψ) = ψ", "основной коммутатор"),
("â<sub>1</sub> = αâ<sub>2</sub> + β*â<sub>2</sub><sup>†</sup>", "α,β из непрерывности ψ и ψ̇ при скачке ω₁→ω₂", "преобразование Боголюбова"),
("n₂ = ⟨0₁∣â<sub>2</sub><sup>†</sup>â<sub>2</sub>∣0₁⟩", "∫∣â<sub>2</sub>ψ<sub>0</sub><sup>(1)</sup>(x)∣²dx = ∣β∣²", "число рождённых квантов"),
("∣α∣² − ∣β∣² = 1", "—", "сохранение коммутатора")]),
]

def build(out):
    S = []
    S.append(Paragraph("Словарь соответствий: бра-кет ↔ интегралы Ландау–Давыдова", T))
    S.append(Paragraph("Перевод дираковской записи на язык волновых функций ψ(x), интегралов и "
        "дифференциальных операторов (Ландау–Лифшиц т. III; Давыдов гл. I–III). Правило: "
        "⟨…∣…⟩ всегда означает интеграл ∫(…)*(…)dx. Разделы 8–9: вывод коэффициента "
        "Боголюбова и формулы Степановского без бра-кет.", Sub))
    for title, rows in SEC:
        S.append(Paragraph(title, H)); S.append(tbl(rows))

    S.append(Paragraph("7. Замечание: КМ частицы и мода поля", H))
    S.append(tbl([
        ("координата x", "амплитуда моды q<sub>k</sub> (ψ(t) в ф.(4) Степановского)", "«координата» осциллятора моды"),
        ("ψ(x)", "Ψ(q<sub>k</sub>)", "волновая функция осциллятора моды"),
        ("ψ(x₁,x₂,…)", "Ψ(q<sub>k1</sub>,q<sub>k2</sub>,…) → функционал Ψ[φ(x)]", "много мод поля"),
        ("∣n⟩ — уровень частицы", "∣n<sub>k</sub>⟩ — n фотонов в моде k", "фотон = квант осциллятора моды")],
        head=("КМ частицы", "Теория поля (Швингер)", "Смысл")))
    S.append(Spacer(1, 4))
    S.append(Paragraph("Таблица раздела 6 применяется к каждой моде отдельно с заменой x → q<sub>k</sub>.", N))

    # ---------- 8. вывод beta ----------
    S.append(Paragraph("8. Коэффициент Боголюбова из гауссова интеграла", H))
    S.append(Paragraph("Осциллятор моды с массой m. До скачка частота ω₁, после — ω₂. Скачок мгновенный, "
        "поэтому волновая функция в момент скачка не успевает измениться: это «внезапное "
        "возмущение» (Ландау–Лифшиц т. III, §41). Система остаётся в старом вакууме", N))
    S.append(Paragraph("ψ<sub>0</sub><sup>(1)</sup>(q) = (mω₁/πħ)<sup>1/4</sup> exp(−mω₁q²/2ħ),", EQ))
    S.append(Paragraph("но теперь его надо разложить по собственным функциям нового гамильтониана "
        "ψ<sub>n</sub><sup>(2)</sup>(q) (гауссиана другой ширины). Вероятность найти n квантов:", N))
    S.append(Paragraph("w<sub>n</sub> = ∣c<sub>n</sub>∣²,   c<sub>n</sub> = ∫ψ<sub>n</sub><sup>(2)</sup>(q)·ψ<sub>0</sub><sup>(1)</sup>(q) dq.", EQ))
    S.append(Paragraph("<b>Шаг 1. Вакуум → вакуум.</b> Интеграл двух гауссиан:", N))
    S.append(Paragraph("c<sub>0</sub> = (ω₁ω₂)<sup>1/4</sup>(m/πħ)<sup>1/2</sup>∫exp(−m(ω₁+ω₂)q²/2ħ)dq "
        "= (ω₁ω₂)<sup>1/4</sup>√(2/(ω₁+ω₂)) = 1/√α,", EQ))
    S.append(Paragraph("где α = (ω₁+ω₂)/(2√(ω₁ω₂)) = ½(√(ω₂/ω₁)+√(ω₁/ω₂)). Значит вероятность "
        "остаться в вакууме w<sub>0</sub> = 1/α &lt; 1 при ω₁ ≠ ω₂.", N))
    S.append(Paragraph("<b>Шаг 2. Нечётные уровни.</b> ψ<sub>0</sub><sup>(1)</sup> чётна по q, "
        "а ψ<sub>2k+1</sub><sup>(2)</sup> нечётны, поэтому c<sub>2k+1</sub> = 0. "
        "Кванты рождаются только парами: n = 0, 2, 4, …", N))
    S.append(Paragraph("<b>Шаг 3. Чётные уровни</b> (интеграл Гаусса–Эрмита):", N))
    S.append(Paragraph("w<sub>2k</sub> = (1/α)·(2k)!/(4<sup>k</sup>(k!)²)·(β/α)<sup>2k</sup>,   "
        "β = ½(√(ω₂/ω₁)−√(ω₁/ω₂)).", EQ))
    S.append(Paragraph("<b>Шаг 4. Среднее число квантов — одним гауссовым интегралом.</b> "
        "Действуем новым понижающим оператором на старый вакуум. Так как "
        "dψ<sub>0</sub><sup>(1)</sup>/dq = −(mω₁/ħ)q·ψ<sub>0</sub><sup>(1)</sup>,", N))
    S.append(Paragraph("â<sub>2</sub>ψ<sub>0</sub><sup>(1)</sup> = √(mω₂/2ħ)·(q + (ħ/mω₂)d/dq)ψ<sub>0</sub><sup>(1)</sup> "
        "= √(mω₂/2ħ)·(1 − ω₁/ω₂)·q·ψ<sub>0</sub><sup>(1)</sup>.", EQ))
    S.append(Paragraph("Подставляем в n₂ = ∫∣â<sub>2</sub>ψ<sub>0</sub><sup>(1)</sup>∣²dq и берём "
        "единственный интеграл ∫q²∣ψ<sub>0</sub><sup>(1)</sup>∣²dq = ħ/(2mω₁):", N))
    S.append(Paragraph("n₂ = (mω₂/2ħ)(1−ω₁/ω₂)²·ħ/(2mω₁) = (ω₂−ω₁)²/(4ω₁ω₂) = β².", EQ))
    S.append(Paragraph("Итог: коэффициент Боголюбова β — это просто несовпадение ширин двух гауссиан. "
        "При ω₁ = ω₂ получаем β = 0. Тождество α² − β² = 1 проверяется подстановкой.", BOX))
    S.append(tbl([
        ("w<sub>0</sub>", "0.98990", "0.98990"),
        ("w<sub>2</sub>", "9.950·10<sup>−3</sup>", "9.950·10<sup>−3</sup>"),
        ("w<sub>4</sub>", "1.500·10<sup>−4</sup>", "1.500·10<sup>−4</sup>"),
        ("w<sub>odd</sub>", "~10<sup>−33</sup> (нуль счёта)", "0"),
        ("Σw<sub>n</sub>", "1.0000000", "1"),
        ("n₂ = Σn·w<sub>n</sub>", "0.0205149", "β² = 0.0205149"),
        ("∫∣â<sub>2</sub>ψ<sub>0</sub><sup>(1)</sup>∣²dq", "0.0205149", "(√ε−1)²/(4√ε) = 0.0205149")],
        head=("Величина (ε = 1.77)", "Численный интеграл", "Формула")))

    # ---------- 9. Степановский ----------
    S.append(Paragraph("9. Степановский (4)–(11) на языке волновых функций", H))
    S.append(Paragraph("Ю.П. Степановский, «Sonoluminescence and black holes as sources of squeezed light» "
        "(ННЦ ХФТИ, памяти Д.В. Волкова). Ниже — тот же вывод без операторов рождения и уничтожения.", N))
    S.append(tbl([
        ("(4)–(5) ψ̈ + ω²(t)ψ = 0,  ω² = c²k²/ε(t)",
         "Мода поля — осциллятор с координатой q (= ψ у Степановского). Квантово: "
         "iħ∂Ψ/∂t = [−(ħ²/2)∂²/∂q² + ½ω²(t)q²]Ψ(q,t)",
         "уравнение Шрёдингера для волновой функции моды"),
        ("(6) скачок ε → 1: e<sup>−iω₁t</sup> → ((√ε+1)/2)e<sup>−iω₂t</sup> + ((√ε−1)/2)e<sup>+iω₂t</sup>",
         "Ψ(q) в момент скачка не меняется (внезапное возмущение, ЛЛ §41); меняется гамильтониан и его собственные функции",
         "сшивка ψ и ψ̇ при скачке; ω₂/ω₁ = √ε"),
        ("(8) разложение поля по a, a<sup>†</sup>",
         "разложение Ψ по собственным функциям ψ<sub>n</sub>(q) осциллятора",
         "у Степановского — операторная форма"),
        ("(9) связь старых и новых операторов (Боголюбов)",
         "c<sub>n</sub> = ∫ψ<sub>n</sub><sup>(2)</sup>ψ<sub>0</sub><sup>(1)</sup>dq — перекрытия гауссиан разной ширины",
         "то же содержание, раздел 8"),
        ("(10) n₂ = ⟨0∣a<sub>2</sub><sup>†</sup>a<sub>2</sub>∣0⟩",
         "n₂ = Σn∣c<sub>n</sub>∣² = ∫∣â<sub>2</sub>ψ<sub>0</sub><sup>(1)</sup>∣²dq",
         "среднее число квантов после скачка"),
        ("(11) n₂ = (√ε−1)²/(4√ε)",
         "(ω₂−ω₁)²/(4ω₁ω₂) при ω₂/ω₁ = √ε",
         "0.0205 на моду при ε = 1.77")],
        head=("Степановский", "Язык Ландау–Давыдова", "Комментарий")))
    S.append(Spacer(1, 5))
    S.append(Paragraph("Физический смысл. После скачка старая гауссиана не является стационарным состоянием "
        "нового осциллятора: её ширина периодически «дышит» с частотой 2ω₂. Это и есть сжатое "
        "состояние (squeezed vacuum), о котором пишет Степановский, а парность w<sub>2k</sub> — "
        "следствие чётности гауссианы. Формулу (7) статьи в таблицу не выношу: её содержание "
        "нужно сверить по оригиналу.", N))

    def foot(c, d):
        c.setFont("DV", 8); c.setFillColor(HexColor("#777777"))
        c.drawCentredString(A4[0]/2, 12*mm, "{0}".format(c.getPageNumber()))
    SimpleDocTemplate(out, pagesize=A4, leftMargin=17*mm, rightMargin=17*mm,
                      topMargin=16*mm, bottomMargin=18*mm,
                      title="Бра-кет ↔ интегралы").build(S, onFirstPage=foot, onLaterPages=foot)
    print("Создан: {0}".format(out))

if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "bra_ket_dictionary.pdf")
