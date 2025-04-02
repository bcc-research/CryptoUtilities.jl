begin 
    for T in [BinaryPoly64, BinaryPoly32, BinaryPoly16, BinaryPoly8]
        a = T(1)
        b = T(2)
        a*b
    end

    for T in [BinaryElem16, BinaryElem32, BinaryElem128]
        a = T(1)
        b = T(2)
        a*b
    end
end