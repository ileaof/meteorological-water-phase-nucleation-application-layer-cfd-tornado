# Auditoria offline do damping superior

Data: 2026-09-10.

A condição superior multiplica `w` quatro vezes por passo. Nas quatro faces do
perfil, os multiplicadores vão de 0,95 na face mais baixa afetada a 1,0 na face
superior. A aplicação não contém `dt`; sua taxa equivalente por segundo depende
do número de passos.

Na janela comum 2790,579551–3300,023494 s:

| medida | 600 m | 300 m | razão 300/600 |
|---|---:|---:|---:|
| passos nativos sobrepostos | 1.015 | 1.283 | 1,264 |
| chamadas fracionárias | 4.057,0 | 5.132,0 | 1,265 |
| taxa nominal na face inferior (s^-1) | 0,4085 | 0,5167 | 1,265 |
| RMS mediano de `w` na camada (m/s) | 0,2443 | 0,2185 | 0,894 |
| máximo de `|w|` na camada (m/s) | 2,7037 | 2,3693 | 0,876 |
| fração mediana do proxy de energia removível por chamada | 0,04234 | 0,04550 | 1,075 |

O operador é nominalmente 26,5% mais forte por segundo em 300 m. Os estados
arquivados também têm menor velocidade vertical na camada e uma distribuição
que perderia 7,5% mais do proxy de energia a cada reaplicação. Essas relações
mostram que a comparação 300/600 não é uma sensibilidade espacial pura.

A perda calculada é contrafactual: aplica-se o operador uma vez a um estado já
arquivado após o passo. Ela não mede a perda causal acumulada, porque tendências
entre chamadas podem repor `w`. O proxy omite densidade e volumes de controle
das faces; somente sua fração dentro de cada caso é comparada.

Artefatos: `outputs/top_damping_audit_20260910/summary.json` e
`snapshot_operator_effect.csv`.
