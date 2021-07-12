## import modules
using Plots
using FFTW
using Images
using HTTP
gr()

## USAF 1951 resolution test chart
usaf= load(HTTP.URI("https://tinyurl.com/c4jtz7cf"))
Ig = Float64.(Gray.(usaf))
ug = sqrt.(Ig)
(M,N) = size(usaf)


L = 0.004	# image size
Δu = L/M
u = -L/2:Δu:(L/2)-Δu

## Functions
function circular_aperture(L,M,fnum)

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

## Interactive
begin

	fnum_slide = @bind fnum_interact Slider(20:1:40, default=30, show_value=true)

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

	width_slide = @bind w_interact Slider(10:0.2:100.0, default=50.0, show_value=true)

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

## Visualize the compoenents

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


plot(p3,p4, layout = (1,2),size = (1200,550), leg = false)
