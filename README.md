# Simulador de bucha oscilante para pá Kaplan — metodologia e referências

## Modelo escolhido

O script representa uma **bucha de deslizamento oscilante no pivô/trunnion de uma pá do rotor Kaplan**. Essa é uma hipótese plausível para estudar atrito e desgaste no cubo; não é uma identificação confirmada do material ou do modelo de mancal das duas turbinas do dataset.

O desgaste usa a forma de Archard:

```text
V = K × W × s / H
h = V / (D × L)
```

em que `V` é volume removido, `K` coeficiente de desgaste, `W` carga normal, `s` distância de deslizamento, `H` dureza e `h` profundidade média estimada. Para uma bucha que oscila entre `-θ` e `+θ`, o código assume distância por ciclo completo de ida e volta `s_ciclo = 4 × r × θ`, com `r = D/2`. O desgaste por passo é calculado em função de uma carga proxy e do número de ciclos assumido no intervalo.

## O que vem do dataset e o que é hipótese

O dataset oferece pressão média/mínima/máxima dos lados de abertura e fechamento, estado de presença de glicerol e, dependendo da turbina, temperatura do óleo hidráulico ou pressão do sistema. Não traz unidades de pressão, geometria/dimensões do mancal, carga normal no pivô, dureza/material, movimento angular das pás, quantidade de ciclos de articulação, medição de desgaste, inspeções ou tempos até falha.

No código, `abs(p_closeAvg - p_openAvg)` é normalizada e usada como **proxy relativo** para variar uma carga nominal assumida. A amplitude intraminuto das pressões é usada como **proxy relativo** de atividade para variar a frequência assumida de articulação. Essas relações não foram identificadas como leis físicas do mancal. A temperatura do óleo hidráulico é preservada no arquivo, mas não entra na equação de desgaste, pois não é a temperatura medida na bucha. O fator de glicerol é neutro (`1.0`) por padrão: o artigo associado relata queda do atrito após injeção, mas não fornece um coeficiente de desgaste comparável.

Os padrões do script — carga de referência, diâmetro, comprimento, dureza, amplitude angular, ciclos/hora e limite de desgaste — são **valores demonstrativos**, e precisam ser substituídos por dados do fabricante, desenho, literatura específica do material ou ensaio. Se o usuário não dispõe de `K` medido em ensaio, o modo `--nominal-life-days` calcula um `K` efetivo para atingir o limite de desgaste na vida nominal escolhida. Essa opção gera cenários para demonstração; ela não descobre a vida real. Se houver ensaio tribológico, use `--k-archard` para fornecer o coeficiente medido.

## Execução

Dependências: Python 3, `numpy` e `pandas`.

```bash
pip install numpy pandas
python main.py \
  --zip kaplan_dataset.zip \
  --out saidas_bucha \
  --turbines Turbine_1 Turbine_2 \
  --stride 15 \
  --scenarios 3 \
  --nominal-life-days 365
```

O ZIP é baixado da fonte oficial se não estiver no caminho indicado. Para resolução minuto a minuto, use `--stride 1`; isso aumenta tempo e tamanho da saída. Para executar em Colab, carregue o script e o ZIP e adapte os caminhos:

```python
!pip install -q numpy pandas
!python /content/main.py \
  --zip /content/kaplan_dataset.zip \
  --out /content/saidas_bucha \
  --stride 15 --scenarios 3 --nominal-life-days 365
```

Os números do exemplo abaixo apenas demonstram a sintaxe. **Não são especificações Kaplan publicadas**:

```bash
python main.py \
  --zip kaplan_dataset.zip --out saidas_exemplo \
  --turbines Turbine_1 --stride 15 --scenarios 3 \
  --nominal-life-days 365 \
  --load-ref-n 10000 --bearing-diameter-mm 150 --bearing-length-mm 100 \
  --hardness-mpa 50 --theta-amplitude-deg 3 --cycles-per-hour 2 \
  --wear-limit-um 500 --glycerol-factor 1.0
```

Alternativamente, se já existir um coeficiente de desgaste obtido de ensaio, use `--k-archard <valor>` no lugar de `--nominal-life-days`. Não use os dois argumentos simultaneamente.

Saídas: `rul_sintetico_bucha_Turbine_1.csv.gz`, `rul_sintetico_bucha_Turbine_2.csv.gz` e `resumo_parametros_e_cenarios.csv`. Os campos `wear_depth_sim_um`, `RUL_sim_hours` e `is_failure_event_sim` são gerados pelo simulador. São rótulos sintéticos condicionados às hipóteses registradas no resumo; não são medições da usina.

## Uso acadêmico recomendado

Use o código como prova de conceito e análise de sensibilidade. Separe trajetórias completas (`scenario_id`) entre treino e teste, em vez de dividir aleatoriamente pontos adjacentes. Compare cenários ao variar carga nominal, número de ciclos, amplitude, limite de desgaste, fator de glicerol e `K`. Não apresente o erro medido contra rótulos sintéticos como desempenho validado em campo. Para um RUL real, seriam necessários ensaios ou inspeções da bucha, critério de fim de vida do fabricante, materiais e dimensões reais e calibração de `K` sob água/glicerol, cargas e oscilação representativas.

## Referências e o que cada uma sustenta

1. **Sandström, L.-J.; Berglund, K.; Marklund, P.; Simmons, G. F. (2025). “Friction Monitoring in Kaplan Turbines.” _Machines_, 13(4), 313.** [Artigo](https://www.mdpi.com/2075-1702/13/4/313) · [DOI](https://doi.org/10.3390/machines13040313). Fonte do contexto do dataset: mancais das pás no cubo, acionamento hidráulico, sinais de pressão e presença de glicerol. O artigo diz que as superfícies expostas podem ser lubrificadas a água ou operar a seco e não publica material/geometria exatos das duas unidades.
2. **Pereira, P.; Schmitt, P.; Riahi, K.; Müller-Brodmann, M. (2009). “Application of self-lubricating bearings in Kaplan runner hubs.” _The International Journal on Hydropower & Dams_, 16(6).** [Página do artigo](https://www.hydropower-dams.com/articles/application-of-self-lubricating-bearings-in-kaplan-runner-hubs/). Sustenta que buchas autolubrificantes são uma opção para cubos Kaplan em arranjos cheios de óleo, cheios de água ou secos. O texto integral pode exigir compra.
3. **Berglin, E. J. et al. (2025). _Self-Lubricating Bushing Testing for John Day Dam Kaplan Turbine Replacement_. PNNL-37577.** [Resumo e relatório](https://www.pnnl.gov/publications/self-lubricating-bushing-testing-john-day-dam-kaplan-turbine-replacement). Caso de outra turbina Kaplan: compara bucha ranhurada de bronze, lubrificada a óleo, com opções autolubrificantes de vários fornecedores; mede atrito e perda de material. Não identifica a unidade do dataset sueco.
4. **Khadke, B. (2023). _Pilot Study: Test Rig for Self-lubricating Kaplan Runner Bearings_.** [Registro e texto](https://www.diva-portal.org/smash/record.jsf?pid=diva2:1737138). Estudo de arranjo de um cubo Kaplan de uma unidade Vattenfall para uma bancada; descreve posições de mancais internos, radial externo e axial/de escora. É um estudo de caso, não uma configuração universal.
5. **U.S. Bureau of Reclamation (2016). _Alignment of Vertical Shaft Hydro Units_, FIST 2-1.** [PDF](https://www.usbr.gov/power/data/fist/fist2_1/fist2-1_9.16.pdf). Referência para distinguir mancais-guia radiais e mancais de escora axiais em conjuntos hidrelétricos verticais e para exemplos de construções desses mancais.
6. **Shutin, D. et al. (2023). “Method for On-Line Remaining Useful Life and Wear Prediction for Adjustable Journal Bearings Utilizing a Combination of Physics-Based and Data-Driven Models: A Numerical Investigation.” _Lubricants_, 11(1), 33.** [Artigo](https://www.mdpi.com/2075-4442/11/1/33). Exemplo de combinação de Archard, lubrificação hidrodinâmica e modelo de dados para mancal de deslizamento; não é uma validação para a turbina deste dataset.
