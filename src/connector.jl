# MaterialStream connector - directed stream (no Kirchhoff flow balance)
@connector function MaterialStream(; name)
    @variables begin
        n(t)
        z(t)
        h(t)
    end
    ODESystem(Equation[], t, [n, z, h], []; name)
end