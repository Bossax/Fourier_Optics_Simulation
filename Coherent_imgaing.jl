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
Simulate an imaging system consisting of a lens with a diameter of 1 inch and a focal length of 600 mm. The LED has 0.5 $\mu$m wave length illuminating an object of USAF 1951 test chart as displayed above. Compute samping criteria ...
"

# ╔═╡ b830713c-65e5-4b1e-9450-48e281d75cee
begin
	flens = 125
	Dlens = 5
	λ = 0.5e-6
	fnum = flens/Dlens
	max_Δu = round(λ*fnum/2, sigdigits = 3)
	max_L = round(M*max_Δu, sigdigits = 3)
	
	@show "max sample interval = $max_Δu m; max side length = $max_L m"
	
end

# ╔═╡ bf3d3bd7-5a3a-4ffe-8ff2-fbb7b6c8a712
md"
Choose side length of the object plane to be 0.01 m
"

# ╔═╡ 5a802c3f-f211-4fc4-9d73-f0641f9d3af7
begin
	L = 0.625e-3# image size
	Δu = L/M
	u = -L/2:Δu:(L/2)-Δu
	
	# # plot a physical object image
	heatmap(u,u,reverse(Ig, dims = 1), xlabel = "m", ylabel = "m",
			color = :grays, aspect_ratio = 1, cbar = true,
		title = "Irridance", titlefontsize = 10, reverse = true, lims = (u[1], u[end]))
	
end

# ╔═╡ 5349202f-7fd2-4915-b386-ae0a014f0dad
md"
Lens MTF
"


# ╔═╡ 79a48318-2364-4ca7-a096-4b65101139d4
begin
	zxp = flens	
	wxp = Dlens/2
	f0 = wxp/(λ*zxp)
	fcoord = -1/(2*Δu):1/L:1/(2*Δu)
	
	# MTF of the exit pupil
	# create meshgrid
	fgrid = [(j,i) for j in fcoord, i in fcoord]
	H = zeros(N,N)
	
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
	
	surface(fcoord./10^5,fcoord./10^5,H, xlabel = "10^5 cycle/m", ylabel = "10^5 cycle/m",zlabel = "Magnitude",color = :grays, size = (500,500), aspect_ratio = 1,
			title = "Magnitude Transfer Function", titlefontsize = 10)
		
end

# ╔═╡ 571b5f6a-8d83-41d6-bddf-4098c40447af
md"
Resulting Image
"

# ╔═╡ 6fb31ace-9784-47cf-aeb9-8a873c2a1468
begin
	U = fftshift(fft(ug))
	U_filtered = U.*H
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
	
	width_slide = @bind w_interact Slider(1:1:20, default=5, show_value=true)
	
	md"""mechanical slit width = $width_slide mm"""

end

# ╔═╡ ce35bdcb-094b-4e35-b574-775aee1cd042
begin
	
	angle_slide = @bind angle_interact Slider(0:1:90, default=90, show_value=true)
	
	md"""mechanical slit angle = $angle_slide degree"""

end

# ╔═╡ 539c43b2-b7ab-4482-862d-c5c8b5f21989
function slit(wide,angle,L,sample,isslit = true)
	du = L/sample
	u = -L/2:du:(L/2)-du
	u_center = -L/2+du:du:L/2-2du
	# create meshgrid of centroids
	grid = [(j,i) for j in u_center, i in u_center]
	
	Slit = zeros(sample,sample)
	if isslit & angle == 90
		line_coord = [(0,k) for k in -sample:sample-1]
		
	else
		return Slit.+1
	end
	
	
	
	
end

# ╔═╡ 86f1308a-4203-4872-858b-3dc52a3c2b9b


# ╔═╡ 10f97269-6b1c-42ba-afaf-da6dd140e004
begin
	# visualize mechanical slit
	#fy = tan(θ)fx
	if is_mech_slit
		slit_width = width_slide
		fcoord;
	else
		heatmap(fcoord,fcoord,ones(length(fcoord),length(fcoord)),
		color = :grays, aspect_ratio = 1, cbar = true,
		title = "Mechanical Slit on Fourier Plane", titlefontsize=10,
			xlabel = "10^5 f(cycle/m)", ylabel =  "10^5 f(cycle/m)", lims = (fcoord[1],fcoord[end]))
		
	end
	
end

# ╔═╡ Cell order:
# ╟─9a879c8d-96c6-4d77-9df8-38530bfa62ce
# ╟─2ef1f8ae-ca62-11eb-01c0-95c668069bb0
# ╟─82e8b964-57b8-4732-985e-d251543090b7
# ╟─d951f366-f8d5-4d2b-b2ea-7f832733f5d1
# ╟─45241e51-d53e-40fa-b306-236ed05901e9
# ╟─62baadb6-9ef8-4c3d-b447-e9622550953a
# ╟─b830713c-65e5-4b1e-9450-48e281d75cee
# ╟─bf3d3bd7-5a3a-4ffe-8ff2-fbb7b6c8a712
# ╠═5a802c3f-f211-4fc4-9d73-f0641f9d3af7
# ╟─5349202f-7fd2-4915-b386-ae0a014f0dad
# ╠═79a48318-2364-4ca7-a096-4b65101139d4
# ╟─571b5f6a-8d83-41d6-bddf-4098c40447af
# ╟─6fb31ace-9784-47cf-aeb9-8a873c2a1468
# ╟─940f8715-826d-4a6d-8567-e03dbb912483
# ╟─b2ce2744-eda6-407f-92ee-8c3e82d797e8
# ╟─696b732a-91d7-4a8f-a9e6-51e889b77d44
# ╟─83692351-63ec-4933-b932-ddde1908e924
# ╟─45b832c1-191a-4254-8394-32183ee9e62c
# ╠═d28475ba-d9aa-4ac3-a031-e7f030190943
# ╟─ce35bdcb-094b-4e35-b574-775aee1cd042
# ╟─539c43b2-b7ab-4482-862d-c5c8b5f21989
# ╠═86f1308a-4203-4872-858b-3dc52a3c2b9b
# ╠═10f97269-6b1c-42ba-afaf-da6dd140e004
