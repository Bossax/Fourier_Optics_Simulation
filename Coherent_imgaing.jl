### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 9a879c8d-96c6-4d77-9df8-38530bfa62ce
# import modules
begin 
	using Plots
	using FFTW
	using Images
	using PlutoUI
	using Plots.PlotMeasures

	gr()
end

# ╔═╡ 2ef1f8ae-ca62-11eb-01c0-95c668069bb0
md"
## Coherent Imaging System

Paraxial f number

$\begin{gather}
f/ \# = \frac{z}{D_{ExP}}
\end{gather}$

Cut-off frequency of a circular aperture

$\begin{gather} f_0 = \frac{1}{2\lambda (f/\#)} \end{gather}$


Nyquist frequency

$\begin{gather}
f_N = \frac{1}{2\Delta u}
\newline \newline
f_o \le \frac{f_N}{2}
\newline
\end{gather}$

Sampling criterion

$\begin{gather}
\Delta u \le \frac{\lambda (f/\#)}{2}
\end{gather}$


### USAF 1951 resolution test chart
---
"

# ╔═╡ d951f366-f8d5-4d2b-b2ea-7f832733f5d1
begin
	usaf= load("usaf_test_chart.jpeg")
	Ig = Float64.(Gray.(usaf));
	ug = sqrt.(Ig);
	usaf
end

# ╔═╡ 45241e51-d53e-40fa-b306-236ed05901e9
begin	
	(res_x,res_y) = size(usaf)
	@show "Image resolution $res_x × $res_y pixels"
 end

# ╔═╡ 62baadb6-9ef8-4c3d-b447-e9622550953a
md"
Simulate an imaging system consisting of a lens with a diameter of 12.5 mm and a focal length of 125 mm. The LED has 0.5 $\mu$m wave length illuminating an object of USAF 1951 test chart as displayed above 
"

# ╔═╡ d469ca7c-e024-4b26-86f7-05b9ddb9ad41
let 
	
	f = 125e-3
	d = 12.5e-3
	λ = 0.5e-6
	f_num = f/d
	Δu = λ*f_num/2

	L = res_x*Δu
	L_dis =round(L*1e3,sigdigits= 4)
	
 	u = -L/2:Δu:L/2 - Δu
	v = u
	
	heatmap(u.*1000,v.*1000,reverse(Ig, dims = 1 ), xlabel = "mm", ylabel = 			"mm",color = :grays, aspect_ratio = 1, cbar = false,
			title = "Ideal Image", titlefontsize = 10,
			lims = (u[1], v[end]))
	
	
	# @show "Sampling period = $(Δu*1e6) μm, Side length = $L_dis mm"
	
end

# ╔═╡ Cell order:
# ╟─9a879c8d-96c6-4d77-9df8-38530bfa62ce
# ╟─2ef1f8ae-ca62-11eb-01c0-95c668069bb0
# ╟─d951f366-f8d5-4d2b-b2ea-7f832733f5d1
# ╟─45241e51-d53e-40fa-b306-236ed05901e9
# ╟─62baadb6-9ef8-4c3d-b447-e9622550953a
# ╠═d469ca7c-e024-4b26-86f7-05b9ddb9ad41
