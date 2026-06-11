# MaterialStream connector - directed stream (no Kirchhoff flow balance)
@connector function MaterialStream(; name)
    @variables begin
        n(t)
        z(t)
    end
    ODESystem(Equation[], t, [n, z], []; name)
end