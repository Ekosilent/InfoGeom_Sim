# InfoGeom-Sim
Interactive simulation of information-induced spacetime geometry and geodesic flows — a 2D toy model.

Real-time SOR Poisson solver with Dirichlet boundary conditions in the browser demo (pure HTML5/JS); a verified DST Dirichlet solver (2nd-order accurate, MMS-tested) powers the Python version and the figures.

> [!IMPORTANT]
> This is a didactic toy model inspired by emergent-gravity ideas (Jacobson, Verlinde, Van Raamsdonk), not a derivation of gravity nor a physical theory. See *Cosa è / cosa non è*.

*(Documentazione tecnica e note in italiano qui sotto. Technical content below is in Italian.)*

---

## Demo

*   **Browser (zero installazione):** apri `docs/index.html`, oppure attiva GitHub Pages (*Settings → Pages → Source: branch (main), cartella (/docs)*) e ottieni un link condivisibile. Trascini le masse, regoli $A$, $B$, $\alpha$, attivi le geodetiche in tempo reale.
*   **Desktop (Python):**
    ```bash
    pip install -r requirements.txt
    python src/interactive_desktop.py
    ```
    *Versione completa:* profilo 1D del potenziale, linee di livello, sorgenti trascinabili, geodetiche on/off. Richiede un backend interattivo (eseguire come script, non in notebook inline).

---

## Il modello

Scalare di informazione e accoppiamento alla geometria (forma concettuale):
$$I(x)=a\,\rho(x)+b|\nabla\rho(x)|^{2}+c\mathcal{F}[\phi(x)], \quad G_{\mu\nu}=\alpha\,I_{\mu\nu}$$

Ciò che il codice implementa davvero è la sola riduzione scalare, statica, in campo debole e 2D, con $\boxed{c=0}$:
$$I=A\rho+B|\nabla\rho|^{2}, \quad \nabla^{2}\Phi=\alpha\,I, \quad h_{00}=-2\Phi, \quad \ddot{r}=-\nabla\Phi$$

L'equazione di Poisson è l'analogo diretto del limite newtoniano ($\nabla^{2}\Phi=4\pi G\rho$), con lo scalare $I$ come sorgente. Le "geodetiche" sono traiettorie di particelle-test newtoniane (non l'equazione geodetica completa).

<img width="1280" height="640" alt="field" src="https://github.com/user-attachments/assets/9783bdde-851e-414a-93db-b8ae8aa2dd2f" />


## Metodo numerico

*   **Poisson con condizioni di Dirichlet** ($\Phi=0$ al bordo) risolto via DST-I, che diagonalizza esattamente il Laplaciano a 5 punti. Niente artefatti periodici, nessun caso singolare ($k=0$).
*   **Variante browser:** SOR red-black con *warm start* per l'interattività in tempo reale.
*   **Traiettorie** integrate con Runge-Kutta su campo di forza ($-\nabla\Phi$) interpolato.

---

## Verifica

Il solver è **verificato** con il **metodo delle soluzioni manufatte (MMS)**: si impone una soluzione esatta che soddisfa le condizioni al contorno (BC) di Dirichlet, si calcola analiticamente il termine sorgente, lo si dà al solver e si confronta.

**Risultato:** convergenza del secondo ordine (pendenza $\approx 1.99$), errore $L^2$ relativo fino a $4\times10^{-6}$.

Questo separa gli errori di modello da quelli di implementazione: il solver è corretto, quindi ogni discrepanza con la fisica è imputabile alle ipotesi, non a bug numerici.

## Cosa è / cosa non è

| È | Non è |
| :--- | :--- |
| Un toy model 2D didattico | Una derivazione della gravità |
| Un solver di Poisson verificato | Una teoria con predizioni falsificabili |
| Una visualizzazione interattiva | Un modello calibrato su scale fisiche |

Nel codice **non esistono**: il campo ausiliario $\phi$, il tensore $I_{\mu\nu}$, il nucleo non locale $\Theta_{\mu\nu}$, l'equazione tensoriale $G_{\mu\nu}=\alpha I_{\mu\nu}$, il vincolo di conservazione $\nabla^{\mu}I_{\mu\nu}=0$, né le geodetiche nulle (raggi di luce). La documentazione tecnica dentro `docs/` contiene la mappa completa *implementato vs aspirazionale*.

---

## Limiti

*   2D, statico, campo debole; nessun limite di campo forte né dinamica della metrica.
*   Accoppiamento $\nabla^{2}\Phi=\alpha I$ postulato, non derivato da un principio d'azione.
*   Unità arbitrarie; $A, B, \alpha$ non calibrati.
*   Nessuna predizione che differisca da Newton/GR.

<img width="2084" height="729" alt="mms" src="https://github.com/user-attachments/assets/3d5415f0-c718-4286-8a60-13040c11702e" />


## Roadmap verso un modello fisico

1. Un'azione $S[\phi, g]$ da cui derivare $I_{\mu\nu}$ e l'accoppiamento, garantendo automaticamente $\nabla^{\mu}I_{\mu\nu}=0$.
2. Dimostrare un limite che recuperi $\nabla^{2}\Phi=4\pi G\rho$, fissando $\alpha$ in termini di $G$.
3. Almeno una predizione confrontabile con i dati (curve di rotazione, lensing, ritardo temporale).

---

## Struttura del repo

```text
InfoGeom-Sim/
├── README.md
├── LICENSE
├── requirements.txt
├── docs/
│   ├── index.html          # demo browser (sorgente GitHub Pages)
│   └── technical_note.pdf  # nota tecnica (logica, metodo, limiti, verifica)
├── figures/
│   ├── field.png           # mappa del campo h00
│   └── mms.png             # grafico di convergenza MMS (mostrato sopra)
└── src/
    ├── interactive_desktop.py  # demo Python completa
    └── make_figures.py         # rigenera le figure della nota
```

---

## Riferimenti

* T. Jacobson, *Thermodynamics of Spacetime: The Einstein Equation of State*, PRL 75, 1260 (1995), [arXiv:gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004)
* E. Verlinde, *On the Origin of Gravity and the Laws of Newton*, JHEP 04 (2011) 029, [arXiv:1001.0785](https://arxiv.org/abs/1001.0785)
* E. Verlinde, *Emergent Gravity and the Dark Universe*, SciPost Phys. 2, 016 (2017), [arXiv:1611.02269](https://arxiv.org/abs/1611.02269)
* M. Van Raamsdonk, *Building up spacetime with quantum entanglement*, GRG 42, 2323 (2010), [arXiv:1005.3035](https://arxiv.org/abs/1005.3035)
* T. Padmanabhan, *Thermodynamical Aspects of Gravity: New insights*, RPP 73, 046901 (2010), [arXiv:0911.5004](https://arxiv.org/abs/0911.5004)
* J. A. Wheeler, *Information, Physics, Quantum: The Search for Links* (1990)

---

## Licenza

MIT — vedi [`LICENSE`](LICENSE).
