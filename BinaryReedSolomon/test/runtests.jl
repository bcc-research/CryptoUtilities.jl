using BinaryReedSolomon
using BinaryFields
using Test

@testset "FFT -> IFFT -> FFT" begin
    n = 10

    for T in [BinaryElem16, BinaryElem32, BinaryElem128]
        twiddles = compute_twiddles(T, n)
        v = rand(T, 2^n)
        v_fft = copy(v)

        fft!(v_fft; twiddles)
        @test v_fft != v

        ifft!(v_fft; twiddles)
        @test v_fft == v
    end
end

@testset "Reed-Solomon systematic" begin
    n, m = 10, 12

    for T in [BinaryElem16, BinaryElem32, BinaryElem128]
        rs = reed_solomon(T, 2^n, 2^m)

        v = rand(T, 2^n)
        output = encode(rs, v)

        @test v == output[1:length(v)]
    end
end

@testset "Sage comparison" begin
    twiddles = BinaryElem128[
        261638842414339399087820898299661203057
        130069497421973758441410450219780457337
        130069497421973758441410450219780457327
        321833370528025984051659201621984161951
        321833370528025984051659201621984161945
        321833370528025984051659201621984161923
        321833370528025984051659201621984161925
        12427004391475801277045897380390817389
        12427004391475801277045897380390817391
        12427004391475801277045897380390817385
        12427004391475801277045897380390817387
        12427004391475801277045897380390817381
        12427004391475801277045897380390817383
        12427004391475801277045897380390817377
        12427004391475801277045897380390817379
    ]

    v = BinaryElem128[48843935073701397021918627474152975110
       257371465678647658219914792930422930533
       197874898248752057839214693713406247745
       86301329031543269357031453671330949739
       245592208151890074913079678553060805151
       191477208903117015546989222243599496680
       92830719409229016308089219817617750833
       264528954340572454088312978462893134650
       158998607558664949362678439274836957424
       187448928532932960560649099299315170550
       177534835847791156274472818404289166039
       307322189246381679156077507151623179879
       117208864575585467966316847685913785498
       332422437295611968587046799211069213610
       109428368893056851194159753059340120844
       197947890894953343492199130314470631788
    ]

    fft!(v; twiddles)

    expected_output = BinaryElem128[
       158767388301301679479875672416174428978
       314045034570696402167150862131636536652
       284497668870731088162348333798389710619
       97193893883131285058688322382264085283
       205661608125885827099961349024782346648
       319854111638988388244315927516461386689
       98163024092465731168779447832503918216
       72461851808861674126157547294435083817
       284672699909608556571358413615868654015
       310357233410493697565822377542976784819
       194488171086938407232562634984109949841
       26083141281753905375688425869148524863
       144700278945341024867563900932218299937
       303726834571845133663217501483978191357
       228881976351733870473775839456225427817
       41896060989421038344777134899638496709
    ]

    @test v == expected_output
end
