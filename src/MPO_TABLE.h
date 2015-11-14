
{
  // if(L < R)
  I(L:f) x I(S:i) -> I(R:f)
  H(L:f) x I(S:i) -> H(R:f)
  I(L:f) x H(S:i) -> H(R:f) / eval V
  C(L:f) x I(S:i) -> C(R:f)
  I(L:f) x C(S:i) -> C(R:f)
  C(L:f) x C(S:i) -> D(R:f)
  D(L:f) x I(S:i) -> D(R:f)
  I(L:f) x D(S:i) -> D(R:f)
  Q(L:f) x I(S:i) -> Q(R:f)
  I(L:f) x Q(S:i) -> Q(R:f) / eval V
  D(L:f) x C(S:i) -> Q(R:f) / eval V
  C(L:f) x D(S:i) -> Q(R:f) / eval V
  Q(L:f) x C(S:i) -> H(R:f)
  C(L:f) x Q(S:i) -> H(R:f)
  D(L:f) x D(S:i) -> H(R:f) / eval V



C1 | C2 |<C3>| C4 | ** | **
C1 | C2 |<C3>| ** | C5 | **
C1 | C2 |<C3>| ** | ** | C6

C1 | C2 | ** |<C4>| C5 | **
C1 | C2 | ** |<C4>| ** | C6
C1 | ** | C3 |<C4>| C5 | **
C1 | ** | C3 |<C4>| ** | C6
** | C2 | C3 |<C4>| C5 | **
** | C2 | C3 |<C4>| ** | C6
C1 | C2 | ** | ** |<C5>| C6

C1 | ** |<C3>| ** | C5 | C6
** | C2 |<C3>| ** | C5 | C6

C1 | ** | ** |<C4>| C5 | C6
** | C2 | ** |<C4>| C5 | C6
** | ** | C3 |<C4>| C5 | C6

C1 | C2 | ** | ** | ** |<C6>| ** | C8
C1 | ** | C3 | ** | ** | C6 | ** | C8
C1 | ** |<C3>| ** | ** | ** | C7 | C8


// CONJ_D case
CREA_CREA x DESA_DESA
CREB_CREB x DESB_DESB
DESA_DESA x CREA_CREA
DESB_DESB x CREB_CREB

// CONJ_D SPIN_CONJ_D and TRANS
CREA_DESA x CREA_DESA
          x CREB_DESB
          x DESA_CREA (i != j)
          x DESB_CREB (i != j)
CREB_DESB x CREA_DESA
          x CREB_DESB
          x DESA_CREA (i != j)
          x DESB_CREB (i != j)
DESA_CREA x CREA_DESA
          x CREB_DESB
          x DESA_CREA (i != j)
          x DESB_CREB (i != j)
DESB_CREB x CREA_DESA
          x CREB_DESB
          x DESA_CREA (i != j)
          x DESB_CREB (i != j)

// CONJ_D and TRANS
CREA_CREB x DESA_DESB
          x DESB_DESA (i != j)
CREB_CREA x DESA_DESB
          x DESB_DESA (i != j)
CREA_DESB x CREB_DESA (i != j)
          x DESA_CREB
CREB_DESA x CREA_DESB
          x DESB_CREA (i != j)
DESA_CREB x CREA_DESB
          x DESB_CREA (i != j)
DESB_CREA x CREB_DESA (i != j)
          x DESA_CREB
DESA_DESB x CREA_CREB
          x CREB_CREA (i != j)
DESB_DESA x CREA_CREB
          x CREB_CREA (i != j)

}
