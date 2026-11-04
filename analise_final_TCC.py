import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from tkinter import Tk, filedialog
from scipy.optimize import least_squares


# ============================================================
# 1. LEITURA DOS ARQUIVOS DO NOVA
# ============================================================

def ler_arquivo_nova(caminho):

    df = pd.read_csv(
        caminho,
        sep=";",
        decimal=","
    )

    col_f = "Frequency (Hz)"
    col_real = "Z' (Ω)"
    col_imag = "-Z'' (Ω)"

    for coluna in [col_f, col_real, col_imag]:
        if coluna not in df.columns:
            raise ValueError(
                f"Coluna '{coluna}' não encontrada."
            )

    f = pd.to_numeric(
        df[col_f],
        errors="coerce"
    ).to_numpy()

    z_real = pd.to_numeric(
        df[col_real],
        errors="coerce"
    ).to_numpy()

    menos_z_imag = pd.to_numeric(
        df[col_imag],
        errors="coerce"
    ).to_numpy()

    Z = z_real - 1j * menos_z_imag

    mascara = (
        np.isfinite(f)
        & np.isfinite(Z.real)
        & np.isfinite(Z.imag)
        & (f > 0)
    )

    f = f[mascara]
    Z = Z[mascara]

    ordem = np.argsort(f)

    return f[ordem], Z[ordem]


# ============================================================
# 2. MODELOS
# ============================================================

def modelo_R(parametros, f):

    R = parametros[0]

    return np.full(
        len(f),
        R,
        dtype=complex
    )


def modelo_RC(parametros, f):

    R, C = parametros

    omega = 2 * np.pi * f

    return 1 / (
        (1 / R)
        + 1j * omega * C
    )


def modelo_CPE(parametros, f):

    R, Q, alpha = parametros

    omega = 2 * np.pi * f

    Y = (
        (1 / R)
        + Q * (1j * omega) ** alpha
    )

    return 1 / Y


# ============================================================
# 3. RESÍDUOS
# ============================================================

def criar_residuo(Z_modelo, Z_exp):

    return np.concatenate([
        Z_modelo.real - Z_exp.real,
        Z_modelo.imag - Z_exp.imag
    ])


def residuos_R(parametros, f, Z_exp):

    return criar_residuo(
        modelo_R(parametros, f),
        Z_exp
    )


def residuos_RC(parametros, f, Z_exp):

    return criar_residuo(
        modelo_RC(parametros, f),
        Z_exp
    )


def residuos_CPE(parametros, f, Z_exp):

    return criar_residuo(
        modelo_CPE(parametros, f),
        Z_exp
    )


# ============================================================
# 4. MÉTRICAS
# ============================================================

def erro_relativo(Z_exp, Z_modelo):

    return (
        100
        * np.linalg.norm(Z_exp - Z_modelo)
        / np.linalg.norm(Z_exp)
    )


def calcular_aicc(residuos, k):

    residuos = np.asarray(residuos)

    n = len(residuos)

    RSS = np.sum(
        residuos ** 2
    )

    if RSS <= 0:
        return -np.inf

    AIC = (
        n * np.log(RSS / n)
        + 2 * k
    )

    if n <= k + 1:
        return np.nan

    return (
        AIC
        + (
            2 * k * (k + 1)
            / (n - k - 1)
        )
    )


# ============================================================
# 5. AJUSTE R
# ============================================================

def ajustar_R(f, Z_exp):

    R0 = np.mean(
        Z_exp.real
    )

    if R0 <= 0:
        R0 = np.mean(
            np.abs(Z_exp.real)
        )

    return least_squares(
        residuos_R,
        x0=[R0],
        bounds=(
            [1e-12],
            [np.inf]
        ),
        args=(f, Z_exp),
        max_nfev=50000
    )


# ============================================================
# 6. AJUSTE R || C
# ============================================================

def ajustar_RC(f, Z_exp):

    R0 = np.max(
        Z_exp.real
    )

    if (
        not np.isfinite(R0)
        or R0 <= 0
    ):
        R0 = np.max(
            np.abs(Z_exp.real)
        )

    indice_pico = np.argmax(
        -Z_exp.imag
    )

    f_pico = f[indice_pico]

    C0 = (
        1
        / (
            2
            * np.pi
            * R0
            * f_pico
        )
    )

    C0 = max(
        C0,
        1e-18
    )

    return least_squares(
        residuos_RC,

        x0=[
            R0,
            C0
        ],

        bounds=(
            [
                1e-12,
                1e-18
            ],
            [
                np.inf,
                np.inf
            ]
        ),

        args=(
            f,
            Z_exp
        ),

        x_scale="jac",

        max_nfev=50000
    )


# ============================================================
# 7. AJUSTE ROBUSTO R || CPE
# ============================================================

def ajustar_CPE(
    f,
    Z_exp,
    R_inicial,
    C_inicial
):

    candidatos = []

    alphas_iniciais = [
        0.50,
        0.70,
        0.85,
        0.95
    ]

    fatores_Q = [
        0.1,
        1.0,
        10.0
    ]

    C_base = max(
        C_inicial,
        1e-18
    )

    for alpha0 in alphas_iniciais:

        for fator_Q in fatores_Q:

            Q0 = max(
                C_base * fator_Q,
                1e-18
            )

            try:

                resultado = least_squares(
                    residuos_CPE,

                    x0=[
                        R_inicial,
                        Q0,
                        alpha0
                    ],

                    bounds=(
                        [
                            1e-12,
                            1e-18,
                            0.10
                        ],
                        [
                            np.inf,
                            np.inf,
                            1.00
                        ]
                    ),

                    args=(
                        f,
                        Z_exp
                    ),

                    x_scale="jac",

                    max_nfev=100000
                )

                RSS = np.sum(
                    resultado.fun ** 2
                )

                if (
                    resultado.success
                    and np.isfinite(RSS)
                ):
                    candidatos.append(
                        (
                            RSS,
                            resultado
                        )
                    )

            except Exception:
                pass

    if not candidatos:

        raise RuntimeError(
            "Nenhum ajuste CPE convergiu."
        )

    candidatos.sort(
        key=lambda x: x[0]
    )

    return (
        candidatos[0][1],
        len(candidatos)
    )


# ============================================================
# 8. CLASSIFICAÇÃO PRELIMINAR
# ============================================================

def classificar_modelo(
    erro_R,
    erro_RC,
    erro_CPE,
    AICc_R,
    AICc_RC,
    AICc_CPE,
    C,
    Q,
    alpha,
    fc_RC_dentro,
    fc_CPE_dentro
):

    # --------------------------------------------------------
    # CPE degenerado
    # --------------------------------------------------------

    cpe_degenerado = (
        alpha <= 0.1001
        or alpha >= 0.9999
        or Q <= 1.01e-18
    )

    # --------------------------------------------------------
    # Espectro muito mal descrito por todos os modelos
    # --------------------------------------------------------

    menor_erro = min(
        erro_R,
        erro_RC,
        erro_CPE
    )

    if menor_erro > 30:

        return (
            "ANÔMALO / REVISAR",
            cpe_degenerado
        )

    # --------------------------------------------------------
    # Compara R com RC
    # --------------------------------------------------------

    delta_R_RC = (
        AICc_R
        - AICc_RC
    )

    # --------------------------------------------------------
    # Compara RC com CPE
    # --------------------------------------------------------

    delta_RC_CPE = (
        AICc_RC
        - AICc_CPE
    )

    # --------------------------------------------------------
    # CPE
    # Só aceita se:
    # - não degenerado
    # - AICc melhor que RC
    # - erro menor
    # - frequência característica observável
    # --------------------------------------------------------

    if (
        not cpe_degenerado
        and delta_RC_CPE > 10
        and erro_CPE < erro_RC
        and fc_CPE_dentro
    ):

        return (
            "R || CPE",
            cpe_degenerado
        )

    # --------------------------------------------------------
    # RC
    # --------------------------------------------------------

    if (
        delta_R_RC > 10
        and erro_RC < erro_R
    ):

        if fc_RC_dentro:

            return (
                "R || C",
                cpe_degenerado
            )

        else:

            return (
                "R || C (fc fora da faixa)",
                cpe_degenerado
            )

    # --------------------------------------------------------
    # Caso predominantemente resistivo
    # --------------------------------------------------------

    return (
        "R",
        cpe_degenerado
    )


# ============================================================
# 9. SELEÇÃO DOS ARQUIVOS
# ============================================================

root = Tk()
root.withdraw()

arquivos = filedialog.askopenfilenames(
    title=(
        "Selecione os arquivos VS01, VS02 e VS03"
    ),

    filetypes=[
        (
            "Arquivos TXT/CSV",
            "*.txt *.csv"
        ),
        (
            "Todos os arquivos",
            "*.*"
        )
    ]
)

root.destroy()

if not arquivos:

    print(
        "Nenhum arquivo selecionado."
    )

    raise SystemExit

# ============================================================
# 10. SELEÇÃO DA PASTA DE SAÍDA
# ============================================================

root = Tk()
root.withdraw()

destino = filedialog.askdirectory(
    title="Selecione a pasta para salvar os resultados"
)

root.destroy()

if not destino:
    print("Nenhuma pasta selecionada. Operação cancelada.")
    raise SystemExit

pasta_saida = Path(destino) / "Resultados_Analise_Final_TCC"
pasta_saida.mkdir(parents=True, exist_ok=True)

pasta_graficos = pasta_saida / "Graficos"
pasta_graficos.mkdir(parents=True, exist_ok=True)


# ============================================================
# 11. PROCESSAMENTO
# ============================================================

resultados = []

print()

print(
    "=" * 80
)

print(
    "ANÁLISE CONSOLIDADA DO TCC"
)

print(
    "Modelos: R | R||C | R||CPE"
)

print(
    "=" * 80
)


for arquivo in arquivos:

    caminho = Path(
        arquivo
    )

    nome = caminho.stem

    print()

    print(
        "-" * 80
    )

    print(
        f"Arquivo: {nome}"
    )

    try:

        # ----------------------------------------------------
        # Dados
        # ----------------------------------------------------

        f, Z_exp = ler_arquivo_nova(
            caminho
        )


        # ====================================================
        # R
        # ====================================================

        ajuste_R = ajustar_R(
            f,
            Z_exp
        )

        R_R = ajuste_R.x[0]

        Z_R = modelo_R(
            ajuste_R.x,
            f
        )

        erro_R = erro_relativo(
            Z_exp,
            Z_R
        )

        AICc_R = calcular_aicc(
            ajuste_R.fun,
            1
        )


        # ====================================================
        # RC
        # ====================================================

        ajuste_RC = ajustar_RC(
            f,
            Z_exp
        )

        R_RC, C_RC = (
            ajuste_RC.x
        )

        Z_RC = modelo_RC(
            ajuste_RC.x,
            f
        )

        erro_RC = erro_relativo(
            Z_exp,
            Z_RC
        )

        AICc_RC = calcular_aicc(
            ajuste_RC.fun,
            2
        )

        fc_RC = (
            1
            / (
                2
                * np.pi
                * R_RC
                * C_RC
            )
        )

        fc_RC_dentro = (
            np.min(f)
            <= fc_RC
            <= np.max(f)
        )


        # ====================================================
        # CPE
        # ====================================================

        (
            ajuste_CPE,
            n_convergencias
        ) = ajustar_CPE(
            f,
            Z_exp,
            R_RC,
            C_RC
        )

        (
            R_CPE,
            Q_CPE,
            alpha
        ) = ajuste_CPE.x

        Z_CPE = modelo_CPE(
            ajuste_CPE.x,
            f
        )

        erro_CPE = erro_relativo(
            Z_exp,
            Z_CPE
        )

        AICc_CPE = calcular_aicc(
            ajuste_CPE.fun,
            3
        )


        # ====================================================
        # FC DO CPE
        # ====================================================

        try:

            fc_CPE = (
                (
                    1
                    / (
                        R_CPE
                        * Q_CPE
                    )
                )
                ** (
                    1 / alpha
                )
                / (
                    2
                    * np.pi
                )
            )

        except (
            OverflowError,
            ZeroDivisionError
        ):

            fc_CPE = np.inf


        fc_CPE_dentro = (
            np.isfinite(fc_CPE)
            and
            np.min(f)
            <= fc_CPE
            <= np.max(f)
        )


        # ====================================================
        # CLASSIFICAÇÃO
        # ====================================================

        (
            modelo_escolhido,
            cpe_degenerado
        ) = classificar_modelo(
            erro_R,
            erro_RC,
            erro_CPE,
            AICc_R,
            AICc_RC,
            AICc_CPE,
            C_RC,
            Q_CPE,
            alpha,
            fc_RC_dentro,
            fc_CPE_dentro
        )


        # ====================================================
        # TERMINAL
        # ====================================================

        print(
            f"R:      erro = {erro_R:.3f}%"
        )

        print(
            f"R||C:   erro = {erro_RC:.3f}%"
        )

        print(
            f"R||CPE: erro = {erro_CPE:.3f}%"
        )

        print(
            f"alpha = {alpha:.4f}"
        )

        print(
            f"CPE degenerado = {cpe_degenerado}"
        )

        print(
            f"Modelo preliminar = {modelo_escolhido}"
        )


        # ====================================================
        # GRÁFICO
        # ====================================================

        plt.figure(
            figsize=(7, 6)
        )

        plt.plot(
            Z_exp.real,
            -Z_exp.imag,
            "o",
            markersize=4,
            label="Experimental"
        )

        plt.plot(
            Z_R.real,
            -Z_R.imag,
            "-",
            label="R"
        )

        plt.plot(
            Z_RC.real,
            -Z_RC.imag,
            "-",
            label="R || C"
        )

        plt.plot(
            Z_CPE.real,
            -Z_CPE.imag,
            "--",
            label="R || CPE"
        )

        plt.xlabel(
            "Z' (Ω)"
        )

        plt.ylabel(
            "-Z'' (Ω)"
        )

        plt.title(
            nome
        )

        plt.legend()

        plt.grid(
            True
        )

        plt.tight_layout()

        plt.savefig(
            pasta_graficos
            / f"{nome}_modelos.png",
            dpi=300
        )

        plt.close()


        # ====================================================
        # DIAGRAMAS DE BODE
        # ====================================================

        # Modelos e dados experimentais
        curvas = [
            ("Experimental", Z_exp, "o"),
            ("R", Z_R, "-"),
            ("R || C", Z_RC, "-"),
            ("R || CPE", Z_CPE, "--")
        ]

        # Bode - módulo da impedância
        plt.figure(figsize=(8, 5))

        for legenda, Z, estilo in curvas:
            plt.loglog(
                f,
                np.abs(Z),
                estilo,
                markersize=4,
                label=legenda
            )

        plt.xlabel("Frequência (Hz)")
        plt.ylabel("|Z| (Ω)")
        plt.title(f"{nome} - Bode (Módulo)")
        plt.legend()
        plt.grid(True, which="both")
        plt.tight_layout()

        plt.savefig(
            pasta_graficos / f"{nome}_Bode_Modulo.png",
            dpi=300
        )

        plt.close()

        # Bode - ângulo de fase
        plt.figure(figsize=(8, 5))

        for legenda, Z, estilo in curvas:
            fase = np.degrees(np.angle(Z))

            plt.semilogx(
                f,
                fase,
                estilo,
                markersize=4,
                label=legenda
            )

        plt.xlabel("Frequência (Hz)")
        plt.ylabel("Ângulo de fase (°)")
        plt.title(f"{nome} - Bode (Fase)")
        plt.legend()
        plt.grid(True, which="both")
        plt.tight_layout()

        plt.savefig(
            pasta_graficos / f"{nome}_Bode_Fase.png",
            dpi=300
        )

        plt.close()

        # ====================================================
        # RESULTADO
        # ====================================================

        resultados.append({

            "Arquivo":
                nome,

            "N pontos":
                len(f),

            "f min (Hz)":
                np.min(f),

            "f max (Hz)":
                np.max(f),

            # R
            "R modelo R (ohm)":
                R_R,

            "Erro R (%)":
                erro_R,

            "AICc R":
                AICc_R,

            # RC
            "R RC (ohm)":
                R_RC,

            "C RC (pF)":
                C_RC * 1e12,

            "fc RC (Hz)":
                fc_RC,

            "fc RC na faixa":
                fc_RC_dentro,

            "Erro RC (%)":
                erro_RC,

            "AICc RC":
                AICc_RC,

            # CPE
            "R CPE (ohm)":
                R_CPE,

            "Q CPE":
                Q_CPE,

            "alpha":
                alpha,

            "fc CPE (Hz)":
                fc_CPE,

            "fc CPE na faixa":
                fc_CPE_dentro,

            "Erro CPE (%)":
                erro_CPE,

            "AICc CPE":
                AICc_CPE,

            "CPE degenerado":
                cpe_degenerado,

            "Convergencias CPE":
                n_convergencias,

            # Comparações
            "Delta AICc R-RC":
                AICc_R - AICc_RC,

            "Delta AICc RC-CPE":
                AICc_RC - AICc_CPE,

            # Resultado
            "Modelo preliminar":
                modelo_escolhido
        })


    except Exception as erro:

        print(
            f"ERRO: {erro}"
        )

        resultados.append({

            "Arquivo":
                nome,

            "Erro":
                str(erro)
        })


# ============================================================
# 12. SALVA RESULTADOS
# ============================================================

df_resultados = pd.DataFrame(
    resultados
)

arquivo_csv = (
    pasta_saida
    / "Resultados_Consolidados_TCC.csv"
)

df_resultados.to_csv(
    arquivo_csv,
    index=False
)


# ============================================================
# 13. RESUMO SIMPLIFICADO
# ============================================================

colunas_resumo = [
    "Arquivo",
    "Erro R (%)",
    "Erro RC (%)",
    "Erro CPE (%)",
    "alpha",
    "Delta AICc R-RC",
    "Delta AICc RC-CPE",
    "CPE degenerado",
    "Modelo preliminar"
]

df_resumo = df_resultados[
    [
        coluna
        for coluna in colunas_resumo
        if coluna in df_resultados.columns
    ]
]

arquivo_resumo = (
    pasta_saida
    / "Resumo_Modelos_TCC.csv"
)

df_resumo.to_csv(
    arquivo_resumo,
    index=False
)


# ============================================================
# 14. FINAL
# ============================================================

print()

print(
    "=" * 80
)

print(
    "ANÁLISE FINALIZADA"
)

print(
    "=" * 80
)

print()

print(
    df_resumo.to_string(
        index=False
    )
)

print()

print(
    "Resultados completos:"
)

print(
    arquivo_csv
)

print()

print(
    "Resumo:"
)

print(
    arquivo_resumo
)

print()

print(
    "Gráficos:"
)

print(
    pasta_graficos
)