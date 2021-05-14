### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 03209d58-df73-4513-8677-46bab1ab424a
# this part is for binder setup
begin
	
	# We set up a new environment for this notebook
	import Pkg
	Pkg.activate(mktempdir())
	
	
	# This is how you add a package:
	Pkg.add("PlutoUI")
	Pkg.add("Plots")
	Pkg.add("FFTW")
	Pkg.add("DSP")
	
end

# ╔═╡ 0e7ca6ce-ae37-11eb-3ec5-099b945b19c4
# import modules
begin 
	using Plots
	using FFTW
	using PlutoUI
	import DSP: unwrap
	using Plots.PlotMeasures

	plotly()
end

# ╔═╡ 655d730b-6dee-45f0-85b1-db04b590a9d8
md"
# Propagation simulation:
Simulate Fresnel and Fraunhofer Diffraction
---
### Frensel diffraction
1. Modelling in frequency domain: Fresnel transfer function

$\begin{gather}
U_2(x,y) = F^{-1}[F[U_1(x,y)H(f_x,f_y)]]\\\\

\textbf{Transer function: } H(f_x,f_y) = e^{jkz} exp[-j\pi \lambda z (f_x^2+f_y^2)]

\end{gather}$

2. Modelling in spatial domain: Fresnel Impulse Response

$\begin{gather}
U_2(x,y) = F^{-1}[F[U_1(x,y)]F[h(x,y)]]\\\\

\textbf{Impulse response: } h(x,y)= \frac{e^{jkz}}{jkz} exp[\frac{jk}{2z}(x^2+y^2)]

\end{gather}$
"

# ╔═╡ 773c805e-65cc-4890-9da6-1108273490d9
md"
Define a propagation transfer function 
"

# ╔═╡ 9bbf46ca-75d0-4d45-8a38-a6635e3a5ce2
function propTF(u1, L, λ, z)
	"""
	Input
	1. source field: u1
	2. sidde length of the plane: L
	3. wavelength: λ
	4. propagation distance: z
	
	Output
	1.field on the observation plane: u2
	
	"""
	j = im
	(M,N) = size(u1)
	dx = L/M 				#sampling distance
	k = 2π/λ 				#wave number
	
	
	f_coord = -1/(2*dx):1/L:1/(2*dx)-1/L 	# frequency coordinate
	H = [exp(-j*π*λ*z*(fx^2 + fy^2)) for fx in f_coord, fy in f_coord]
	H = fftshift(H)
	U1 = ifft(ifftshift(u1))
	
	return ifftshift(ifft(H.*U1))
	
	
end

# ╔═╡ ae2c8d6a-c89f-4ba2-94de-a175c0f48d1a
md"
Define a propagation impulse response
"

# ╔═╡ 0c6f7289-b243-4f97-bdb8-b60d85a16eff
function propIR(u1, L, λ, z)
	"""
	Input
	1. source field: u1
	2. sidde length of the plane: L
	3. wavelength: λ
	4. propagation distance: z
	
	Output
	1.field on the observation plane: u2
	
	"""
	j = im
	(M,N) = size(u1)
	dx = L/M 				#sampling distance
	k = 2π/λ 				#wave number
		
	x_coord = -L/2:dx:L/2-dx 	# spatial coordinate
	h = [1/(j*λ*z)*exp(j*k/(2*z)*(x^2+y^2)) for x in x_coord, y in x_coord]
	H = fftshift(fft(h)).*dx^2
	U1 = fftshift(fft(u1))
	
		
	return ifftshift(ifft(ifftshift(H.*U1)))
	
end

# ╔═╡ 38008805-19a6-4150-b231-9c17aa5f5e3a
md"
### Create a square beam to test the functions

A square source field with dimensions 0.5 m x 0.5 m
"

# ╔═╡ bd6d5054-e2c0-45a9-aecf-60c51790569f
# square funtion
function square(Ls,N,w_x,w_y,A_u = 1.0)
	dxs = Ls/N
	xs = -Ls/2:dxs:Ls/2-dxs
	
	aperture = zeros(N,N)
	
	start_index_x = findlast(xs .<= -w_x)
	stop_index_x = findfirst(xs .>= w_x)
	start_index_y = findlast(xs .<= -w_y)
	stop_index_y = findfirst(xs .>= w_y)
	
	for x in start_index_x:stop_index_x
		for y in start_index_y:stop_index_y
			aperture[y,x] = A_u
		end
	end
	
	return aperture
	
	
end

# ╔═╡ c008a685-9d95-4a15-9430-e69d0cd9d6c9
# circle funtion
function circle(Ls,N,r,A_u = 1.0)
	dxs = Ls/N
	
	# centroids of the pixels
	x_center = -Ls/2+dxs:dxs:Ls/2-2dxs
	# create meshgrid
	grid = [(j,i) for j in x_center, i in x_center]
	aperture = zeros(N,N)
	
	keep_coord = []
	#find (x,y) that have their centroids fall within the radius of the circular aperture
	for x in 1:size(grid,2)
		for y in 1:size(grid,1)
			dis = sqrt(grid[x,y][1]^2+grid[x,y][2]^2)
			dis <= r ? push!(keep_coord,(x,y)) : nothing
		end
		
	end
	
	# assign valued to the selected coordinates
	for k in keep_coord 
		(i,j) = k
		aperture[j,i] = A_u
	end
		
		
	return aperture
	
	
end

# ╔═╡ 15ddd06a-cf27-4333-a06b-ec1c7b17c197
# source parameters
begin 
	Lc = 0.5 			#side length
	Nc = 500 			# samples
	r = 0.051
	
	λ = 0.5e-6
	k =2*π/λ
	z = 2000
	
	u3 = circle(Lc,Nc,r);
	I3 = abs.(u3.^2);

	dxc = Lc/Nc
	x_coordc = -Lc/2:dxc:Lc/2-dxc
	heatmap(x_coordc,x_coordc,u3, xlabel = "m", ylabel = "m",
			color = :grays, size = (300,300), aspect_ratio = 1,
			title = "Circle beam", titlefontsize = 10)
end

# ╔═╡ e6f8522f-a665-4b81-8c43-f91f2807549e
# source parameters
begin 
	Ls = 0.5 			#side length
	N = 250 			# samples
	w_x = 0.051 		# half-width of the aperture
	w_y = 0.051
	
	λc = 0.5e-6
	kc =2*π/λ
	zc = 2000
	
	u1 = square(Ls,N,w_x,w_y);
	I1 = abs.(u1.^2);

	dx = Ls/N
	x_coord = -Ls/2:dx:Ls/2-dx
	heatmap(x_coord,x_coord,u1, xlabel = "m", ylabel = "m",
			color = :grays, size = (300,300), aspect_ratio = 1,
			title = "Square beam", titlefontsize = 10)
end

# ╔═╡ 2a995bf5-125c-44d9-ac24-26b40065d4ca
md"
**Result from the transfer function**
"

# ╔═╡ 814e7c4c-24a1-4c4a-92e5-4dcb1d9ddc97
let
	u2_TF = propTF(u1, Ls, λ, z)
	
	# irridance
	plt1 = heatmap(x_coord,x_coord,abs.(u2_TF.^2), xlabel = "m", ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = false,
			title = "Irridance", titlefontsize = 10,
			lims = (x_coord[1], x_coord[end]))
	# phase
	plt2 = plot(x_coord, unwrap(angle.(u2_TF[convert(Int64,N/2+1),:])),
				xlabel = "m", ylabel = "rad",
				title = "Phase", titlefontsize = 10)
	
	plot(plt1,plt2, layout = (1,2),	size = (500,300), leg = false)
end

# ╔═╡ 67e7a36f-2e64-4efa-9c81-c2908d1f5b28
let
	u2_TF = propTF(u3, Lc, λ, z)
	
	# irridance
	plt1 = heatmap(x_coordc,x_coordc,abs.(u2_TF.^2), xlabel = "m", ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = false,
			title = "Irridance", titlefontsize = 10,
			lims = (x_coordc[1], x_coordc[end]))
	# phase
	plt2 = plot(x_coordc, unwrap(angle.(u2_TF[convert(Int64,Nc/2+1),:])),
				xlabel = "m", ylabel = "rad",
				title = "Phase", titlefontsize = 10)
	
	plot(plt1,plt2, layout = (1,2),	size = (500,300), leg = false)
end

# ╔═╡ 795aad6e-1f7d-4a15-93d5-8325d82b4658
md"
**Result from impulse response**
"

# ╔═╡ d371a0a7-c9d0-459b-bcfe-8e83cc400d43
let
	u2_IR = propIR(u1, Ls, λ, z)
	
	# irridance
	plt1 = heatmap(x_coord,x_coord,abs.(u2_IR.^2), xlabel = "m", ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = false,
			title = "Irridance", titlefontsize = 10,
			lims = (x_coord[1], x_coord[end]))
	
	#phase
	plt2 = plot(x_coord, unwrap(angle.(u2_IR[convert(Int64,N/2+1),:])),
				xlabel = "m", ylabel = "rad",
				title = "Phase", titlefontsize = 10)
	
	plot(plt1,plt2, layout = (1,2),	size = (500,300), leg = false)
end

# ╔═╡ eaf4589e-76f0-4b31-9fd6-188560fd7e3f
let
	u2_IR = propIR(u3, Lc, λ, z)
	
	# irridance
	plt1 = heatmap(x_coordc,x_coordc,abs.(u2_IR.^2), xlabel = "m", ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = false,
			title = "Irridance", titlefontsize = 10,
			lims = (x_coordc[1], x_coordc[end]))
	# phase
	plt2 = plot(x_coordc, unwrap(angle.(u2_IR[convert(Int64,Nc/2+1),:])),
				xlabel = "m", ylabel = "rad",
				title = "Phase", titlefontsize = 10)
	
	plot(plt1,plt2, layout = (1,2),	size = (500,300), leg = false)
end

# ╔═╡ e6f7007d-bbd5-42d2-a3c9-3b60acd83d33
md"
---
## Effect of sampling interval, propagation distance, and aperture size
"

# ╔═╡ f56f2cf9-2828-4ec1-b7d7-f27fbe041df2
begin
	z_slider = @bind z_interact Slider(1000:1000:20000, default=2000, show_value=true)
	md""" Propagation distance = $(z_slider) m"""
end


# ╔═╡ 9e90c16c-5142-42f9-87ac-3958a2ce1bfb
begin
	N_slider = @bind N_interact Slider(100:50:500, default=250, show_value=true)
	md""" Sample = $(N_slider) m"""
end

# ╔═╡ b0357d32-6669-4f2f-aa61-1a81c9612922
begin
	w_slider = @bind w_interact Slider(0.011:0.005:0.081, default=0.051, show_value=true)
	md""" Half-width = $(w_slider) m"""
end

# ╔═╡ 4770d444-db59-4e85-baa0-0d54834afbcd
let
	F_number = round(w_interact^2/(λ*z_interact), sigdigits=2)
	Markdown.MD(Markdown.Admonition("Fresnel Number", "Fresnel Number",
					[md" `` \frac{w^2}{\lambda z}`` =  $F_number , where ``w`` is half-width of the square"]))
	
end

# ╔═╡ d4af9fc6-e198-43bc-8847-5c2fe7b2a590
let
	u1_sq_interact = square(Ls,N_interact,w_interact,w_interact)
	u2_TF = propTF(u1_sq_interact, Ls, λ, z_interact)
	u2_IR = propIR(u1_sq_interact, Ls, λ, z_interact)
	
	dx_sq = Ls/N_interact
	dist_coord = -Ls/2:dx_sq:Ls/2-dx_sq
	
	I_TF = abs.(u2_TF.^2)
	I_IR = abs.(u2_IR.^2)
	
	# irridance
	plt1 = heatmap(dist_coord,dist_coord,I_TF, ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = false,
			title = "Irridance", titlefontsize = 10, 
			lims = (dist_coord[1], dist_coord[end]))
	# x section
	plt2 = plot(dist_coord, I_TF[convert(Int64,N_interact/2+1),:]./ maximum(I_TF[convert(Int64,N_interact/2+1),:]),
				 xlabel = "m",ylabel = "I", aspect_ratio = .5, 
				title = "Magnitude cross-section", titlefontsize = 10)
	
	
	
	
	
	# irridance
	plt3 = heatmap(dist_coord,dist_coord,abs.(u2_IR.^2), ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = false,
			title = "Irridance", titlefontsize = 10,
			lims = (dist_coord[1], dist_coord[end]))
	
	#x section
	plt4 = plot(dist_coord, I_IR[convert(Int64,N_interact/2+1),:]./ maximum(I_IR[convert(Int64,N_interact/2+1),:]),
				xlabel = "m", ylabel = "I",aspect_ratio = .5,
				title = "Magnitude cross-section", titlefontsize = 10)
	
	plot(plt1,plt2,plt3,plt4, layout = @layout([plt1 plt2; plt3 plt4]),
				size = (600,600), leg = false,top_margin = 20px)
end

# ╔═╡ d0aa8b36-1332-47c0-8be9-d33f58395d22
md"
---
### Fraunhofer diffraction

Fresnel Number <<1

$\begin{gather}
U_2(x,y) = \frac{exp(jkz)}{j\lambda z} exp[j\frac{k}{2z}(x_2^2 + y_2^2)] \times \iint U_1(x_1,y_1)exp[-j\frac{2\pi}{\lambda z}(x_2x_1 + y_2y_1)]
\end{gather}$

The side lengths of the source plane and the observation plane are not generally the same.

`` L_2 = \frac{\lambda z}{\Delta x} ``, and `` \Delta x_2 = \frac{\lambda z}{L_1}``

"

# ╔═╡ 950176a0-b190-4908-9561-a7f1cc1949b0
function propFF(u1, L1, λ, z)
	
	"""
	Input
	1. source field: u1
	2. side length of the source plane: L1
	3. wavelength: λ
	4. propagation distance: z
	
	Output
	1.field on the observation plane: u2
	2. side length of the observation plane
	
	"""
	j = im
	(M,N) = size(u1)
	dx1= L1/M 				#sampling distance
	k = 2π/λ 				#wave number
	
	L2 = λ*z/dx1
	dx2 = λ*z/L1
	
	x_coord2 = -L2/2:dx2:L2/2-dx2	# spatial coordinate
	c = [1/(j*λ*z)*exp(j*k/(2*z)*(x^2 +y^2)) for x in x_coord2, y in x_coord2]
	
	return  (c.*ifftshift(fft(u1))*dx1^2, L2)
	
	
end

# ╔═╡ f7c30403-7716-42a1-b7a9-68884456a249
let
	u1FF = square(Ls,N,0.011,0.011);
	
	(u2,L2)= propFF(u1FF,Ls, λc, zc)
	
	dx2 = λc*zc/Ls
	x_coordFF = -L2/2:dx2:L2/2-dx2
	# irridance
	plt1 = heatmap(x_coordFF,x_coordFF,abs.(u2), xlabel = "m", ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = false,
			title = "Field", titlefontsize = 10,
			lims = (x_coordFF[1], x_coordFF[end]))
	# magnitude x section
	plt2 = plot(x_coordFF, abs.(u2[convert(Int64,N/2+1),:]),
				xlabel = "m", ylabel = "irridance",
				title = "Magnitude cross-section", titlefontsize = 10)
	
	plot(plt1,plt2, layout = (1,2),	size = (500,300), leg = false)
	
end

# ╔═╡ 103ad2be-8eb6-400c-93bc-0f73a09c2df7


# ╔═╡ 0ce9ed22-5104-4778-917f-48e8485f9bda
md"
---
"

# ╔═╡ 87d02650-1355-4b3c-b7c7-2a618810f815
md"
## Notes on FFT
"

# ╔═╡ f24801ac-bbc4-4548-b1c2-8515ab22f1ff
begin
	sig = hcat(zeros(1,10),ones(1,20))'
	sig = vcat(sig,zeros(10,1))
	
	p1 = plot(sig, line = :stem, marker = :circle, c = :red, title = "unshifted", xlabel = "n", ylabel = "f(x)")
	p2 = plot(fftshift(sig), line = :stem, marker = :circle, c = :blue, title = "shifted", xlabel = "n", ylabel = "f(x)")

	p = plot(p1,p2,layout = (2,1), size = (600,500), leg = false)
end

# ╔═╡ aedae1ce-87a9-4d26-b232-c392c81c8814
begin
	S1 = fft(sig)
	r1 = plot(0:length(S1)-1,abs.(S1), line = :stem, marker = :circle, c = :red, title = "Magnitude spectrum: unshifted")
	r2 = plot(-length(S1)/2:length(S1)/2-1,abs.(fftshift(S1)), line = :stem, marker = :circle, c = :blue, title = "Magnitude spectrum: shifted")
	plot(r1,r2,layout=(2,1), size = (500,700), leg = false)
	
end

# ╔═╡ afae35d2-3960-432d-b6fc-f2b321fd1dd4
begin
	L_s1 = length(S1)
	ind = convert(Int64,L_s1/2-1)
	S1_folded = vcat(S1[1] ,S1[2:ind]*2)
	s = ifft(ifftshift(fftshift(S1)))		# unshifting
	plot(real(s),line = :stem,marker = :circle, size = (600,300), leg = false)
	
end

# ╔═╡ 7693fcc0-927e-479e-a7f0-b93028ad26f4
length(s)

# ╔═╡ Cell order:
# ╟─655d730b-6dee-45f0-85b1-db04b590a9d8
# ╠═03209d58-df73-4513-8677-46bab1ab424a
# ╠═0e7ca6ce-ae37-11eb-3ec5-099b945b19c4
# ╟─773c805e-65cc-4890-9da6-1108273490d9
# ╠═9bbf46ca-75d0-4d45-8a38-a6635e3a5ce2
# ╟─ae2c8d6a-c89f-4ba2-94de-a175c0f48d1a
# ╠═0c6f7289-b243-4f97-bdb8-b60d85a16eff
# ╟─38008805-19a6-4150-b231-9c17aa5f5e3a
# ╟─bd6d5054-e2c0-45a9-aecf-60c51790569f
# ╟─c008a685-9d95-4a15-9430-e69d0cd9d6c9
# ╟─e6f8522f-a665-4b81-8c43-f91f2807549e
# ╟─15ddd06a-cf27-4333-a06b-ec1c7b17c197
# ╟─2a995bf5-125c-44d9-ac24-26b40065d4ca
# ╟─814e7c4c-24a1-4c4a-92e5-4dcb1d9ddc97
# ╟─67e7a36f-2e64-4efa-9c81-c2908d1f5b28
# ╟─795aad6e-1f7d-4a15-93d5-8325d82b4658
# ╟─d371a0a7-c9d0-459b-bcfe-8e83cc400d43
# ╟─eaf4589e-76f0-4b31-9fd6-188560fd7e3f
# ╟─e6f7007d-bbd5-42d2-a3c9-3b60acd83d33
# ╟─f56f2cf9-2828-4ec1-b7d7-f27fbe041df2
# ╟─9e90c16c-5142-42f9-87ac-3958a2ce1bfb
# ╟─b0357d32-6669-4f2f-aa61-1a81c9612922
# ╟─4770d444-db59-4e85-baa0-0d54834afbcd
# ╟─d4af9fc6-e198-43bc-8847-5c2fe7b2a590
# ╟─d0aa8b36-1332-47c0-8be9-d33f58395d22
# ╠═950176a0-b190-4908-9561-a7f1cc1949b0
# ╠═f7c30403-7716-42a1-b7a9-68884456a249
# ╠═103ad2be-8eb6-400c-93bc-0f73a09c2df7
# ╟─0ce9ed22-5104-4778-917f-48e8485f9bda
# ╟─87d02650-1355-4b3c-b7c7-2a618810f815
# ╟─f24801ac-bbc4-4548-b1c2-8515ab22f1ff
# ╠═aedae1ce-87a9-4d26-b232-c392c81c8814
# ╠═afae35d2-3960-432d-b6fc-f2b321fd1dd4
# ╠═7693fcc0-927e-479e-a7f0-b93028ad26f4
