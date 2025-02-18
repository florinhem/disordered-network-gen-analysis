
import LinearAlgebra
import Statistics
import DataFrames
import PlotlyJS

function energy(r,a_List,theta_eq_degree,beta)
    @assert length(r)===3           #rx,ry,rz
    @assert length(a_List)===4      #a0,a1,a2,a3
    
    #display("1:")
    E=0
    for i in 1:size(a_List,1)
        @assert length(a_List[i])===3
        dE=(LinearAlgebra.dot(r-a_List[i],r-a_List[i])/(3/16)-1)^2
    
        E+=dE
        
    end
   
    theta_eq_radian=theta_eq_degree/360.0*2*pi

    for i in 1:size(a_List,1)
        for j in 1:size(a_List,1)
            if i > j
                dE=2*beta*(LinearAlgebra.dot(r-a_List[i],r-a_List[j])/(3/16)-cos(theta_eq_radian))^2   #Factor 2 because of i>j symmetry
                E+=dE
            end
        end
    end

    display("E, $E")
    


    return E
end

function perturbation_theory(;
    r_List,
    a_List,
    theta_eq_degree_List,
    beta
    )

    energy_theta_List=[]
    for theta_eq_degree in theta_eq_degree_List
        display("theta: $theta_eq_degree")
        energy_r_List=[]
        for r in r_List
            display("r: $r")
            E=energy(r,a_List,theta_eq_degree,beta)
            append!(energy_r_List,E)
        end
        energy_r_mean=Statistics.mean(energy_r_List)
        append!(energy_theta_List,energy_r_mean)
    end

    return energy_theta_List
end

function plot_perturbation_theory(
    theta_eq_degree_List,
    energy_theta_List)

    df_plot = DataFrames.DataFrame(
        x=theta_eq_degree_List,
        y=energy_theta_List,
    )

    PT=PlotlyJS.scatter(
        df_plot,
        x=:x, 
        y=:y, 
        #mode="lines",
        line = PlotlyJS.attr(color = "green")
        ) 

    Delta2_layoutb1b2 = PlotlyJS.Layout(
        title=
            "PT",
        xaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=round.(theta_eq_degree_List,digits=3),
            tickvals=theta_eq_degree_List
        ),
        #=
        yaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=round.(energy_theta_List,digits=3),
            tickvals=energy_theta_List
        ),
        =#
        xaxis_title="theta_eq_degree",
        yaxis_title="energy",
        autosize=false
        
    )

    Plot_PT=PlotlyJS.plot(PT,Delta2_layoutb1b2)

    plot_total_path=raw".\simulations\analysis_plot\PT_11.png"

    PlotlyJS.savefig(Plot_PT,plot_total_path)

    display("plotted")

end

display("------------")
r=[1/4,1/4,1/4]
display(r)
r_List=0.00*rand(3, 1) .+ r
display(r_List)
r_List= [r_List[:,i] for i in 1:size(r_List,2)]
display(r_List)
A=2/sqrt(3)
theta_eq_degree_List=0:10:180
display(theta_eq_degree_List)
#[0,20,40,60,80,90,95,100,105,109.47122063449,115.0,120,125,130,140,160,180]


energy_theta_List=perturbation_theory(
    r_List=r_List,#[r],     #Make random number array
    #a_List=[A*[-1/2,-1/2,-1/2], A*[1/2,1/2,-1/2], A*[1/2,-1/2,1/2], A*[-1/2,1/2,1/2]],      
    #a_List=[[0,0,0], [0,0,1], [2*sqrt(2)/3,0,-1/3], [sqrt(23)/6 , 1/2, -1/3],[sqrt(23)/6 , -1/2, -1/3] ],
    a_List=[[1/2,1/2,0], [1/2,0,1/2], [0,1/2,1/2], [0,0,0]],  
    theta_eq_degree_List=theta_eq_degree_List,
    beta=0.285
    )

plot_perturbation_theory(
    theta_eq_degree_List,
    energy_theta_List
)

display("------------")