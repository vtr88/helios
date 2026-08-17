# helios

Sol e Lua no terminal. Tanabi, SP, é o local padrão.

```sh
make
./helios
./helios 2026-12-21
./helios --ache "Ubatuba, SP"
./helios --local=-23.43389,-45.07111,America/Sao_Paulo
```

Sem data, usa o instante e o dia atuais do local. `--ache` requer internet e
`curl`.

ISC. Cálculos adaptados do SunCalc (BSD-2-Clause); veja
`THIRD_PARTY_LICENSES`.
