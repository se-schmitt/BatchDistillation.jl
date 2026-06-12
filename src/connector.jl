# MaterialStream connector - directed stream (no Kirchhoff flow balance)
# n: molar flow, z: composition vector (length N), h: molar enthalpy
@connector function MaterialStream(; name, N)
    @variables begin
        n(t)
        z(t)[1:N]
        h(t)
    end
    ODESystem(Equation[], t, [n, z, h], []; name)
end