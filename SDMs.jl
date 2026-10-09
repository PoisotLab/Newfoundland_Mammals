using SpeciesDistributionToolkit
const SDT = SpeciesDistributionToolkit
using CairoMakie
using Statistics
using PrettyTables

# Import data

    #Polygon for SDM region
    esri_polygon = getpolygon(PolygonData(ESRI, Places))
    amerinord = esri_polygon["Continent"=>"North America"]

        #BoundingBox
        extent = SDT.boundingbox(amerinord)

    #Download 19 bioclim data from WorldClim
    world_bioclim = RasterData(WorldClim2, BioClim)
    layers(world_bioclim)[1:5]
    
        #Get variables' text description
        [
            layerdescriptions(world_bioclim)[l]
            for l in layers(world_bioclim)
        ]

# Get 19 layers for all variables in the bbox
layers_bioclim = Vector{SDMLayer{Float32}}(world_bioclim, extent)

    #Visualisation
    f = Figure(; size = (600, 300))
    ax = Axis(f[1, 1]; aspect = DataAspect())
    heatmap!(ax, layers[1]; colormap = :Greys)
    hidespines!(ax)
    hidedecorations!(ax)
    f
mask!(layers, amerinord)
