using Plots

# I am just using my hands to input the values because why not
keepProbs = [0.01,0.05,0.1]

f_tr  = []
f_ang = 

err_tr = []
err_54 = []
err_65 = []

open("waveletErrors_transverse.txt","r") do f
    line = 0
    while ! eof(f)
        s = readline(f)
        if length(s)>3
            if s[1:3] == "ERR"
                global err_tr
                error = parse(Float64, s[8:end])
                err_tr = vcat(err_tr,[error])
            end
        end
        line += 1
    end
end

open("waveletErrors_54.txt","r") do f
    line = 0
    while ! eof(f)
        s = readline(f)
        if length(s)>3
            if s[1:3] == "ERR"
                global err_54
                error = parse(Float64, s[8:end])
                err_54 = vcat(err_54,[error])
            end
        end
        line += 1
    end
end

open("waveletErrors_65.txt","r") do f
    line = 0
    while ! eof(f)
        s = readline(f)
        if length(s)>3
            if s[1:3] == "ERR"
                global err_65
                error = parse(Float64, s[8:end])
                err_65 = vcat(err_65,[error])
            end
        end
        line += 1
    end
end


err_tr_a1 = reshape(err_tr[1:75],(3,25))';    err_tr_a2 = reshape(err_tr[76:150],(3,25))'; 
err_tr_a3 = reshape(err_tr[151:225],(3,25))';    err_tr_a4 = reshape(err_tr[226:300],(3,25))'; 

err_54_a1 = reshape(err_54[1:75],(3,25))';    err_54_a2 = reshape(err_54[76:150],(3,25))'; 
err_54_a3 = reshape(err_54[151:225],(3,25))';    err_54_a4 = reshape(err_54[226:300],(3,25))'; 

err_65_a1 = reshape(err_65[1:75],(3,25))';    err_65_a2 = reshape(err_65[76:150],(3,25))'; 
err_65_a3 = reshape(err_65[151:225],(3,25))';    err_65_a4 = reshape(err_65[226:300],(3,25))'; 

println("MAX ERROR TRA2 P05: "*string(maximum(err_tr_a2[:,1])))
println("MAX ERROR TRA4 P05: "*string(maximum(err_tr_a4[:,1])))

println("MAX ERROR TRA2 P10: "*string(maximum(err_tr_a2[:,1])))
println("MAX ERROR TRA4 P10: "*string(maximum(err_tr_a4[:,1])))

println("MAX ERROR 54A2 P05: "*string(maximum(err_54_a2[:,1])))
println("MAX ERROR 54A4 P05: "*string(maximum(err_54_a4[:,1])))

println("MAX ERROR 54A2 P10: "*string(maximum(err_54_a2[:,1])))
println("MAX ERROR 54A4 P10: "*string(maximum(err_54_a4[:,1])))

println("MAX ERROR 65A2 P05: "*string(maximum(err_65_a2[:,1])))
println("MAX ERROR 65A4 P05: "*string(maximum(err_65_a4[:,1])))

println("MAX ERROR 65A2 P10: "*string(maximum(err_65_a2[:,1])))
println("MAX ERROR 65A4 P10: "*string(maximum(err_65_a4[:,1])))
#l = @layout  [a b c]
#p3 = plot(1:25,hcat(err_tr_a2[:,2],err_tr_a2[:,3],err_tr_a4[:,2],err_tr_a4[:,3]))
#p2 = plot(1:25,hcat(err_54_a2[:,2],err_54_a2[:,3],err_54_a4[:,2],err_54_a4[:,3]))
#p1 = plot(1:25,hcat(err_65_a2[:,2],err_65_a2[:,3],err_65_a4[:,2],err_65_a4[:,3]))
#ylabel!("% Error")

#plot(p1,p2,p3,layout=l)
#ylims!((0,0.5))
#xlabel!("Frequency Index")
