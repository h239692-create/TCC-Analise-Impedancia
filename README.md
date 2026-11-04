# Análise de Espectroscopia de Impedância — TCC

**Modelagem e Automação da Extração de Parâmetros Elétricos de Sensores Químicos Baseados em Semicondutores Orgânicos**

**Autor:** Henrique Mariano Freitas  
**Instituição:** Universidade Estadual de Campinas (UNICAMP)  
**Unidade:** Faculdade de Engenharia Elétrica e de Computação (FEEC)  
**Linguagem:** Python

## 1. Descrição do projeto

Este repositório contém a rotina computacional desenvolvida como parte do Trabalho de Conclusão de Curso (TCC) em Engenharia Elétrica na Universidade Estadual de Campinas (UNICAMP).

O programa tem como objetivo automatizar a análise de dados experimentais obtidos por espectroscopia de impedância eletroquímica (*Electrochemical Impedance Spectroscopy* — EIS), aplicada à caracterização elétrica de sensores químicos baseados em semicondutores orgânicos.

A rotina realiza a leitura dos dados experimentais, o ajuste de modelos de circuitos equivalentes, a extração de parâmetros elétricos, a avaliação quantitativa dos ajustes e a geração automática de gráficos e tabelas.

## 2. Modelos de circuitos equivalentes

São considerados três modelos:

### Modelo resistivo (R)

Representa um comportamento predominantemente resistivo:

$$
Z(\omega)=R
$$

### Modelo resistor em paralelo com capacitor (R ∥ C)

Representa um comportamento resistivo-capacitivo ideal:

$$
Z(\omega)=\frac{1}{\frac{1}{R}+j\omega C}
$$

### Modelo resistor em paralelo com elemento de fase constante (R ∥ CPE)

Representa um comportamento resistivo-capacitivo não ideal:

$$
Z(\omega)=\frac{1}{\frac{1}{R}+Q(j\omega)^\alpha}
$$

Em que:

- \(R\): resistência elétrica (Ω);
- \(C\): capacitância elétrica (F);
- \(Q\): coeficiente do elemento de fase constante (CPE);
- \(\alpha\): expoente adimensional do CPE;
- \(\omega=2\pi f\): frequência angular (rad/s);
- \(j\): unidade imaginária.

## 3. Metodologia computacional

O ajuste dos modelos é realizado por mínimos quadrados não lineares, utilizando a função `least_squares` da biblioteca SciPy.

O vetor de resíduos é construído a partir das diferenças entre as componentes real e imaginária das impedâncias experimentais e calculadas.

A avaliação dos modelos considera:

- Erro relativo da impedância complexa;
- Critério de Informação de Akaike corrigido (AICc);
- Comparação entre modelos com diferentes quantidades de parâmetros;
- Verificação da frequência característica em relação à faixa experimental;
- Identificação de soluções degeneradas do CPE;
- Classificação preliminar do circuito equivalente mais adequado.

O ajuste do modelo R ∥ CPE utiliza múltiplas estimativas iniciais, buscando reduzir a dependência da solução em relação aos valores de partida.

## 4. Requisitos

Para executar o programa, é necessário ter o Python instalado, juntamente com as bibliotecas:

- NumPy
- Pandas
- Matplotlib
- SciPy

O programa também utiliza o módulo `tkinter` para a seleção interativa dos arquivos experimentais e da pasta de destino.

As dependências externas estão relacionadas no arquivo `requirements.txt`.

## 5. Instalação e execução

### 5.1. Obter o código

Na página do repositório, clique em **Code → Download ZIP** e extraia os arquivos em uma pasta de sua preferência.

### 5.2. Instalar as dependências

Abra um terminal na pasta do projeto e execute:

```bash
python -m pip install -r requirements.txt
```

Recomenda-se utilizar um ambiente virtual Python para instalar as dependências.

### 5.3. Executar o programa

Execute:

```bash
python analise_final_TCC.py
```

O programa abrirá uma janela para selecionar os arquivos experimentais e, em seguida, outra janela para escolher a pasta onde os resultados serão armazenados.

## 6. Formato dos dados de entrada

O programa foi desenvolvido para analisar arquivos `.txt` ou `.csv` exportados pelo software NOVA, contendo dados de espectroscopia de impedância.

Os arquivos devem apresentar as seguintes colunas:

| Coluna | Descrição |
|---|---|
| `Frequency (Hz)` | Frequência em hertz |
| `Z' (Ω)` | Componente real da impedância |
| `-Z'' (Ω)` | Componente imaginária negativa da impedância |

A rotina espera arquivos com separador de campos `;` e separador decimal `,`.

Os dados são organizados em ordem crescente de frequência antes do processamento.

## 7. Resultados gerados

Após o processamento, o programa cria uma pasta denominada `Resultados_Analise_Final_TCC` no local selecionado pelo usuário.

A estrutura de saída é:

```text
Resultados_Analise_Final_TCC/
├── Resultados_Consolidados_TCC.csv
├── Resumo_Modelos_TCC.csv
└── Graficos/
    ├── nome_do_arquivo_modelos.png
    ├── nome_do_arquivo_Bode_Modulo.png
    └── nome_do_arquivo_Bode_Fase.png
```

### Resultados_Consolidados_TCC.csv

Contém os parâmetros ajustados dos modelos, os erros relativos, os valores de AICc, as frequências características e as informações utilizadas na classificação dos espectros.

### Resumo_Modelos_TCC.csv

Apresenta uma síntese dos indicadores de ajuste e do modelo preliminarmente selecionado para cada espectro.

### Gráficos

Para cada arquivo experimental processado com sucesso, são gerados:

- **Diagrama de Nyquist:** comparação das componentes real e imaginária da impedância;
- **Diagrama de Bode — módulo:** módulo da impedância em função da frequência;
- **Diagrama de Bode — fase:** ângulo de fase da impedância em função da frequência.

Os gráficos apresentam os dados experimentais e as respostas calculadas para os três modelos.

## 8. Aplicação acadêmica

A rotina foi desenvolvida para apoiar a análise e a interpretação de espectros de impedância de sensores químicos orgânicos, permitindo comparar diferentes representações elétricas e automatizar a extração de parâmetros.

A classificação fornecida pelo programa deve ser interpretada em conjunto com o comportamento físico dos sensores e com as características dos dados experimentais.

## 9. Autor e instituição

**Henrique Mariano Freitas**  
Graduação em Engenharia Elétrica  
Faculdade de Engenharia Elétrica e de Computação — FEEC  
Universidade Estadual de Campinas — UNICAMP  
Campinas, São Paulo, Brasil — 2026

Projeto desenvolvido no âmbito do Trabalho de Conclusão de Curso.

## 10. Observações

O código é disponibilizado para consulta e reprodutibilidade acadêmica.

Os resultados dependem da qualidade dos dados experimentais, da faixa de frequências analisada e da adequação dos modelos elétricos considerados.

A disponibilização pública do repositório não estabelece, por si só, uma licença de reutilização ou redistribuição do código.
