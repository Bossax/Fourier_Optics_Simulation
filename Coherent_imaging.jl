### A Pluto.jl notebook ###
# v0.14.8

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

# ╔═╡ 290554ae-64b7-4588-8014-1c2e635fa78b
# this part is for binder setup

begin
	
	# We set up a new environment for this notebook
	import Pkg
	Pkg.activate(mktempdir())
	
	
	# This is how you add a package:
	Pkg.add("PlutoUI")
	Pkg.add("Plots")
	Pkg.add("FFTW")
	Pkg.add("Images")
	
	
end


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
# Coherent Imaging System:


An imaging system consists of lenses to scale object and relay it onto a plane. This operation results in magnification and diffraction which is nature of converging lenses. Since the converging lens applies phase modulation to the wavefront, with the Fresnel propagation law, the result of the imaging system involves Fraunhofer diffraction pattern, which is a Fourier Transform of the source field, esentially pupil functions if other transmittance functions from  phase modulations are excluded. 

*Key concepts*

1. The ideal image produced by a diffraction-limited optical system is a scaled and inverted version of the object (geometrical optics)
2. The effect of diffraction is to convolve that ideal image with the Fraunhofer pattern of the lens pupil.


$\begin{gather}
U_i(x,y) = h(x,y) * U_g(x,y)
\newline
U_g(x,y) = \frac{1}{M} U_o(\xi/M, \eta/M)
\end{gather}$

 and ``h`` is the impulse response of the exit pupil. The magnitude transfer funciton, hence


$\begin{gather}
H(f_X,f_Y) = F[h] = P(-\lambda zf_x, \lambda zf_y)
\end{gather}$

### Sampling in coherent imaging simulation
Paraxial f number

$\begin{gather}
f/ \# = \frac{z}{D_{ExP}}
\end{gather}$

Cut-off frequency of a circular aperture (ideal image)

$\begin{gather} f_0 = \frac{1}{2\lambda (f/\#)} \end{gather}$


Nyquist sampling requency (imposed by the sampled image imported into the workspace)

$\begin{gather}
f_N = \frac{1}{2\Delta u}
\end{gather}$

To capture the entire bandwidth of the ideal image,

$\begin{gather}

f_o \le \frac{f_N}{2}
\newline
\end{gather}$

Sampling criterion

$\begin{gather}
\Delta u \le \frac{\lambda (f/\#)}{2}
\end{gather}$

The size of the image on the image plane is dependent on the sample interval and the number of pixels of the object 

"

# ╔═╡ 82e8b964-57b8-4732-985e-d251543090b7
md"
### USAF 1951 resolution test chart
---
"

# ╔═╡ d951f366-f8d5-4d2b-b2ea-7f832733f5d1
begin
	usaf= load("usaf_test_chart.jpeg")
	Ig = Float64.(Gray.(usaf))
	ug = sqrt.(Ig)
	@show "Loading sampled image"
	
end

# ╔═╡ 45241e51-d53e-40fa-b306-236ed05901e9
begin	
	(M,N) = size(usaf)
	@show "Image resolution $M × $N pixels"
 end

# ╔═╡ 62baadb6-9ef8-4c3d-b447-e9622550953a
md"
Simulate an imaging system consisting of a lens with a diameter of 1 inch and a focal length of 125 mm. The LED has 0.5 $\mu$m wave length illuminating an object of USAF 1951 test chart as displayed above. Compute samping criteria ...
"

# ╔═╡ b830713c-65e5-4b1e-9450-48e281d75cee
begin
	flens = 125 # mm
	Dlens = 5 # mm
	λ = 0.5e-6
	fnum = flens/Dlens
	max_Δu = round(λ*fnum/2, sigdigits = 3)
	max_L = round(M*max_Δu, sigdigits = 3)
	
	@show "max sample interval = $max_Δu m; max side length = $max_L m"
	
end

# ╔═╡ bf3d3bd7-5a3a-4ffe-8ff2-fbb7b6c8a712
md"
Choose side length of the object plane to be 0.000625 m
"

# ╔═╡ 5a802c3f-f211-4fc4-9d73-f0641f9d3af7
begin
	L = 0.625e-3# image size
	Δu = L/M
	u = -L/2:Δu:(L/2)-Δu
	
	# # plot a physical object image
	usaf
	
end

# ╔═╡ 5349202f-7fd2-4915-b386-ae0a014f0dad
begin

function circular_aperture(L,M,fnum)
	# exit pupil
	# create meshgrid
	δu = L/M
	f0 = (λ*fnum)^-1
	fcoord = -1/(2*δu):1/L:1/(2*δu)
	fgrid = [(j,i) for j in fcoord, i in fcoord]
	H = zeros(M,M)
	
	keep_coord = []
	
	for fx in 1:size(fgrid,2)
		for fy in 1:size(fgrid,1)
			fdis = sqrt(fgrid[fx,fy][1]^2+fgrid[fx,fy][2]^2)
			fdis <= f0 ? push!(keep_coord,(fx,fy)) : nothing
		end
		
	end
	# assign valued to the selected coordinates
	for k in keep_coord 
		(i,j) = k
		H[j,i] = 1
	end
	return H
end
end

# ╔═╡ 79a48318-2364-4ca7-a096-4b65101139d4
begin
	
	circular_lens1 = circular_aperture(L,M,fnum)
	fcoord = -1/(2*Δu):1/L:1/(2*Δu)
	surface(fcoord./10^5,fcoord./10^5,circular_lens1, xlabel = "10^5 cycle/m", ylabel = "10^5 cycle/m",zlabel = "Magnitude",color = :grays, size = (500,500), aspect_ratio = 1,
			title = "Magnitude Transfer Function", titlefontsize = 10)
		
end

# ╔═╡ 571b5f6a-8d83-41d6-bddf-4098c40447af
md"
Resulting Image
"

# ╔═╡ 6fb31ace-9784-47cf-aeb9-8a873c2a1468
begin
	U = fftshift(fft(ug))
	U_filtered = U.*circular_lens1
	ui = ifft(U_filtered)
# Image Amplitude Spectrum
	p1 = heatmap(fcoord./(10^5), fcoord./(10^5), log.(abs.(U)), 
			aspect_ratio = 1, 
			color = :oslo, colorbar =true, colorbar_scales = :log10, 
			colorbar_title = "Log(magnitude)",
			title = "Image Fourier Amplitude Spectrum", titlefontsize=10,
			xlabel = "10^5 f(cycle/m)", ylabel =  "10^5 f(cycle/m)",
		lims = (fcoord[1]./(10^5), fcoord[end]/(10^5)));
	
	p2= heatmap(fcoord./(10^5), fcoord./(10^5), log.(abs.(U_filtered)), 
		aspect_ratio = 1,
		color = :oslo, colorbar =true, colorbar_scales = :log10,
		colorbar_title = "Log(magnitude)",
		title = "Filtered Amplitude Spectrum", titlefontsize=10,
		xlabel = "10^5 f(cycle/m)", ylabel =  "10^5 f(cycle/m)",
		lims = (fcoord[1]./(10^5), fcoord[end]/(10^5)));
	@show "Computing...."
	
end

# ╔═╡ 940f8715-826d-4a6d-8567-e03dbb912483
plot(p1,p2, layout = (1,2),size = (800,450), leg = false)

# ╔═╡ b2ce2744-eda6-407f-92ee-8c3e82d797e8
heatmap(u,u,reverse(abs.(ui).^2, dims = 1), xlabel = "m", ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = true,
		title = "Blurred Irridance f/# = $fnum", titlefontsize = 10, reverse = true, lims = (u[1], u[end]))

# ╔═╡ 696b732a-91d7-4a8f-a9e6-51e889b77d44
md"
**Optical components**
1. Lens 
2. Mechanical slit
3. Central obscuration
"

# ╔═╡ 83692351-63ec-4933-b932-ddde1908e924
begin
	
	fnum_slide = @bind fnum_interact Slider(5:1:20, default=10, show_value=true)
	
	md""" Lens f/ $fnum_slide """

end

# ╔═╡ 45b832c1-191a-4254-8394-32183ee9e62c
begin
	mech_checkbox = @bind is_mech_slit CheckBox(default=false)
	
	md"""
	Mechanical Slit $mech_checkbox 
	"""
end


# ╔═╡ d28475ba-d9aa-4ac3-a031-e7f030190943
begin
	
	width_slide = @bind w_interact Slider(0.8:0.2:10.0, default=5.0, show_value=true)
	
	md"""mechanical slit width = $width_slide μm"""

end

# ╔═╡ ce35bdcb-094b-4e35-b574-775aee1cd042
begin
	
	angle_slide = @bind angle_interact Slider(0:1:90, default=90, show_value=true)
	
	md"""mechanical slit angle = $angle_slide degree"""

end

# ╔═╡ d64dcbe5-d78a-4bde-801c-2ce9a41adfee
begin
	cobs_checkbox = @bind is_cobs CheckBox(default=false)
	
	md"""
	Central obscuration $cobs_checkbox 
	"""
end

# ╔═╡ 2641cc17-a78b-4b11-b044-a511b8f6f332
begin
	
	r_slide = @bind r_interact Slider(0.4:0.2:2, default=1, show_value=true)
	
	md"""Central obscuration radius % = $r_slide %"""

end

# ╔═╡ db6ecd92-e0a7-4d8d-8822-cda7dcb11026
begin
	r_cobs = L/2*r_interact/100*1000
	@show "Radius = $r_cobs mm"
end


# ╔═╡ 539c43b2-b7ab-4482-862d-c5c8b5f21989
function slit(width,angle,L,sample,isslit = true)
	
	if ~isslit
		return ones((sample,sample))
	end
	
	
	# Domain
	w = width/1e6
	M = sample
	Δx = L/M
	x = -L/2:Δx:L/2-Δx
	slit = zeros((M,M))
	
	# Line
	mid_point = convert(Int64,M/2)
	
	if angle == 90

		x_bound = convert(Int64,floor(w/2/Δx))
		for k in mid_point-x_bound:mid_point+x_bound
			slit[:,k] .= 1
		end
		return slit
	elseif angle == 0
		
		y_bound = convert(Int64,floor(w/2/Δx))
		for k in mid_point-y_bound:mid_point+y_bound
			slit[k,:] .= 1
		end
		return slit
		
	end
	
	θ = angle
	m = round(tan(θ/180*pi),sigdigits = 3)
	C = M/2*(1-m)
	cells =[]
	
	# indices
	for x_now in 1:M
		y = convert(Int64,floor(m*x_now+C))
		if y >= 2
			dist = [(y+g)-m*x_now for g in -1:1]
			(~,k) = findmin(dist)
			push!(cells,(x_now,y))
		else
			dist = [(y+g)-m*x_now for g in 0:1]
			(~,k) = findmin(dist)
			push!(cells,(x_now,y))
		end

	end

	# search radius
	opening_cells = []
	cell_cov = convert(Int,ceil(w/2/Δx))
	perpen_m = -1/m
	x_dis = convert(Int,ceil(w/2*sin(θ/180*pi)/Δx))
	for cell in cells
		local i,j   
		(i,j) = cell
		intercept = j+i/m

		x_lower = i-x_dis
		x_lower < 1 ? x_lower = 1 : nothing
		x_lower > M ? x_lower = M : nothing
		x_upper = i+x_dis
		x_upper < 1 ? x_upper = 1 : nothing
		x_upper > M ? x_upper = M : nothing

		y_lower = convert(Int,floor(perpen_m*x_upper + intercept) )
		y_lower < 1 ? y_lower = 1 : nothing
		y_lower > M ? y_lower = M : nothing
		y_upper = convert(Int,ceil(perpen_m*x_lower + intercept))
		y_upper < 1 ? y_upper = 1 : nothing
		y_upper > M ? y_upper = M : nothing

		# get the indices
		(x_upper, y_lower) ∉ opening_cells ? push!(opening_cells,(x_upper, y_lower)) : nothing
		(x_lower, y_upper) ∉ opening_cells ? push!(opening_cells,(x_lower, y_upper)) : nothing
		#check if it hits the upper edge
		if y_lower == M
			break;
		end


	end
	## fill in the gap

	x_coord =[cell[1] for cell in opening_cells]
	y_coord =[cell[2] for cell in opening_cells]
	max_x = maximum(x_coord)

	p_lower_y = 0
	p_upper_y = M+1

	for x_now in 1:max_x
		cell_set = findall(h -> h == x_now, x_coord )
		y_points = [y_coord[g] for g in cell_set]
		upper_y = maximum(y_points)
		lower_y = minimum(y_points)
		y_dist = upper_y-lower_y

		# case 1 the shaded area is well bounded
		if y_dist != 0
			for k in lower_y:upper_y
				slit[k,x_now] = 1
			end
			p_upper_y = upper_y
		    p_lower_y = lower_y

		else 
			# case 2 no lower bound
			if p_lower_y == 1
				for k in 1:upper_y
					slit[k,x_now] = 1
				end

			# # case 3 no upper bound
			elseif p_upper_y == M
				for k in lower_y:M
					slit[k,x_now] = 1
				end
			end

		end



	end
	
	return slit
end

# ╔═╡ 3e5cfcac-eba3-4595-9df0-4916d24ba4b7
function central_obsuration(L,M,flens,r_ratio, is_obs = false)
	if ~is_obs 
		return ones(M,M)
	end
	r = L/2*r_ratio
	δu = L/M
	u = -L/2:δu:L/2-δu
	grids = [(j,i) for j in u, i in u]
	cobs = ones(M,M)
	
	keep_coord = []
	
	for x in 1:size(grids,2)
		for y in 1:size(grids,1)
			xdist = sqrt(grids[x,y][1]^2+grids[x,y][2]^2)
			xdist <= r ? push!(keep_coord,(x,y)) : nothing
		end
		
	end
	# assign valued to the selected coordinates
	for k in keep_coord 
		(i,j) = k
		cobs[j,i] = 0
	end
	return cobs
end

# ╔═╡ 10f97269-6b1c-42ba-afaf-da6dd140e004
begin
	slit1 = slit(w_interact,angle_interact,L,M, is_mech_slit)
	cobs1 = central_obsuration(L,M,flens,r_interact/100, is_cobs);
	# visualize mechanical slit 
	p5 = heatmap(u.*1000,u.*1000,slit1, aspect_ratio = 1, color = :grays, xlabel = "mm", ylabel = "mm",
		lims = (u[1]*1000, u[M]*1000),size = (300,250))
	p6 = heatmap(u.*1000,u.*1000,cobs1, aspect_ratio = 1, color = :grays, xlabel = "mm", ylabel = "mm",
		lims = (u[1]*1000, u[M]*1000),size = (300,250))
	
	plot(p5,p6, layout = (1,2),size = (1200,550), leg = false)
	
end

# ╔═╡ 150ddf48-6a2d-4671-aaaa-f585eaa0a6ab
begin
	circular_lens_interact = circular_aperture(L,M,fnum_interact)
	U_slit_filtered = U.*slit1.*circular_lens_interact.*cobs1
	ui2 = ifft(U_slit_filtered)
	# Image Amplitude Spectrum
	p3 = heatmap(fcoord./(10^5), fcoord./(10^5), log.(abs.(U_slit_filtered)), 
		aspect_ratio = 1,
		color = :oslo, colorbar =true, colorbar_scales = :log10,
		colorbar_title = "Log(magnitude)",
		title = "Filtered Amplitude Spectrum", titlefontsize=10,
		xlabel = "10^5 f(cycle/m)", ylabel =  "10^5 f(cycle/m)",
		lims = (fcoord[1]./(10^5), fcoord[end]/(10^5)));
	
	p4 = heatmap(u.*1000,u.*1000,reverse(abs.(ui2).^2, dims = 1), xlabel = "mm",
		ylabel = "mm", color = :grays, aspect_ratio = 1, cbar = true,
		title = "Slit Filtered f/# = $fnum_interact", titlefontsize = 10,
		reverse = true, lims = (u[1]*1000, u[end]*1000))
	# resulting image
	@show "Computing Result...."
	
end

# ╔═╡ 189414d4-22d6-42d4-a737-fa9648385388
plot(p3,p4, layout = (1,2),size = (1200,550), leg = false)

# ╔═╡ Cell order:
# ╟─290554ae-64b7-4588-8014-1c2e635fa78b
# ╟─9a879c8d-96c6-4d77-9df8-38530bfa62ce
# ╟─2ef1f8ae-ca62-11eb-01c0-95c668069bb0
# ╟─82e8b964-57b8-4732-985e-d251543090b7
# ╟─d951f366-f8d5-4d2b-b2ea-7f832733f5d1
# ╟─45241e51-d53e-40fa-b306-236ed05901e9
# ╟─62baadb6-9ef8-4c3d-b447-e9622550953a
# ╟─b830713c-65e5-4b1e-9450-48e281d75cee
# ╟─bf3d3bd7-5a3a-4ffe-8ff2-fbb7b6c8a712
# ╟─5a802c3f-f211-4fc4-9d73-f0641f9d3af7
# ╟─5349202f-7fd2-4915-b386-ae0a014f0dad
# ╟─79a48318-2364-4ca7-a096-4b65101139d4
# ╟─571b5f6a-8d83-41d6-bddf-4098c40447af
# ╟─6fb31ace-9784-47cf-aeb9-8a873c2a1468
# ╟─940f8715-826d-4a6d-8567-e03dbb912483
# ╟─b2ce2744-eda6-407f-92ee-8c3e82d797e8
# ╟─696b732a-91d7-4a8f-a9e6-51e889b77d44
# ╟─83692351-63ec-4933-b932-ddde1908e924
# ╟─45b832c1-191a-4254-8394-32183ee9e62c
# ╟─d28475ba-d9aa-4ac3-a031-e7f030190943
# ╟─ce35bdcb-094b-4e35-b574-775aee1cd042
# ╟─d64dcbe5-d78a-4bde-801c-2ce9a41adfee
# ╟─2641cc17-a78b-4b11-b044-a511b8f6f332
# ╟─db6ecd92-e0a7-4d8d-8822-cda7dcb11026
# ╟─539c43b2-b7ab-4482-862d-c5c8b5f21989
# ╟─3e5cfcac-eba3-4595-9df0-4916d24ba4b7
# ╟─10f97269-6b1c-42ba-afaf-da6dd140e004
# ╟─150ddf48-6a2d-4671-aaaa-f585eaa0a6ab
# ╠═189414d4-22d6-42d4-a737-fa9648385388
