using SpeciesDistributionToolkit
const SDT = SpeciesDistributionToolkit
using CairoMakie
using Dates
using Statistics
using PrettyTables

# Import data

    #Polygon for SDM region
    amerinord = getpolygon(PolygonData(ESRI, Places))["Continent"=>"North America"]["Land type" => "Primary land"]

        #BoundingBox
        #extent = SDT.boundingbox(amerinord) (on prend pas ce bbox prcq il est vrm trop gros)
        extent = (left = -168.75, right=13.24, bottom=-8., top=84.4)

    #Download 19 bioclim data from WorldClim
    world_bioclim = RasterData(WorldClim2, BioClim)
    layers(world_bioclim)[1:5]
    
        #Get variables' text description
        [
            layerdescriptions(world_bioclim)[l]
            for l in layers(world_bioclim)
        ]

    # Get 19 layers for all variables in the bbox
    layers_bioclim = Vector{SDMLayer{Float32}}(world_bioclim, amerinord);
    mask!(layers_bioclim, amerinord)

        #Visualisation
        f = Figure(; size = (600, 300))
        ax = Axis(f[1, 1]; aspect = DataAspect())
        heatmap!(ax, layers_bioclim[1]; colormap = :Greys)
        hidespines!(ax)
        hidedecorations!(ax)
        f

# Microtus pennsylvanicus
    # Import GBIF data 
        records1 = GBIF.download("10.15468/dl.z4yb6u") #all datasets exept iNaturalist
        test_records1 = GBIF.download("10.15468/dl.h8vnkf") #iNaturalist dataset
        campagnol = taxon(first(entity(records1)))

        #Get the Phylopic silouhette
        silhouette_campagnol = Phylopic.imagesof(campagnol; items = 1)
        Phylopic.attribution(silhouette_campagnol)

        #Visualisation
        f = Figure()
        ax = Axis(f[1, 1]; aspect = DataAspect())
        poly!(ax, amerinord; color = :grey95)
        scatter!(ax, records1; color = :forestgreen, label = "Occurences M. pennsylvanicus")
        scatter!(
            ax,
            test_records1;
            color = (:purple, 0.4),
            strokecolor = :purple,
            strokewidth = 1,
            marker = :rect,
            markersize = 9,
            label = "Test",
        )
        hidedecorations!(ax)
        hidespines!(ax)
        lines!(ax, amerinord; color = :grey10)
        axislegend(ax; position = :rt, framevisible = false)
        silhouetteplot!(
            ax,
            extent.left - 0.3,
            extent.top + 0.2,
            silhouette_campagnol;
            markersize = 70,
            color = :black,
            )
        f
#Generate pseudo-absences

    #Enlever les occurrences qui n'ont pas de coordonnées spatiales
records1 = filter(e -> !ismissing(place(e)), elements(records1))
records1 = mask(Occurrences(records1), amerinord)


    #Layer with presence data
    presencelayer = SimpleSDMLayers.mask(first(layers_bioclim), Occurrences(records1))


    #Pseudo-absence mask
background = sqrt.(pseudoabsencemask(DistanceToEvent, presencelayer))

    #Sample the pseudo absences (** à modifier biologiquement pour déterminer les paramètres de dispersion du campagnol)
    bgpoints = backgroundpoints(
    nodata(background, x -> !(sqrt(10.0) <= x <= sqrt(100.0))), #pseudo absences crées au minimum à 4 km d'une occurrence et à 30 km d'une occurrence
    2sum(presencelayer), #On génère 2 fois plus de pseudo absences que d'occurrences (à modifier si nécessaire). 2 fois c'est standard
)

    #Visualisation
    f = Figure(; size = (600, 300))
ax = Axis(f[1, 1]; aspect = DataAspect())
poly!(ax, amerinord; color = :grey90, strokecolor = :black, strokewidth = 1)
scatter!(ax, presencelayer; color = :black)
scatter!(ax, bgpoints; color = :red, markersize = 2)
hidedecorations!(ax)
hidespines!(ax)
f

#Setting up the model
model_microtus = SDM(PCATransform, Logistic, layers_bioclim, presencelayer, bgpoints)

#Tweak the hyperparameter (à faire plus tard)

#Variable selection -cross validation
folds = kfold(model_microtus; k = 5)
variables!(model_microtus, ForwardSelection, folds; included = [1]) #adding variables 1 at a time
istrained(model_microtus) #Check that the model is trained
variables(model_microtus) #Inspect the list of trained variables
    #Ce sont les variables 1,3,16 et 18 qui expliquent le mieux le modèle

#Performance evaluation of the model
cval, ctrn = crossvalidate(model_microtus)
cval[1]

    #Tableau de Performance
    ms = [mcc, tpr, fpr, fnr, tnr, ppv, npv]
M = permutedims([
    m(c) for m in ms, c in [cval, ctrn]
]);
M = hcat(["Training", "Validation"], M);
pretty_table(
    M;
    alignment = [:l, :c, :c, :c, :c, :c, :c, :c],
    backend = :markdown,
    column_labels = ["Dataset"; uppercase.(string.(ms))],
    formatters = [fmt__printf("%5.3f", [2, 3, 4, 5, 6, 7, 8])],
)

#Prediction on the testing data
Xₜ = layers_bioclim[test_records1]
yₜ = predict(model_microtus, Xₜ)
accuracy(yₜ, ones(Bool, length(yₜ))) #testing the accuracy

#Understanding the model (à faire ultérieurement)

#Map
P = predict(model_microtus, layers_bioclim; threshold = false)

f = Figure(; size = (600, 300))
ax = Axis(f[1, 1]; aspect = DataAspect())
hm = heatmap!(ax, P; colormap = Reverse(:batlowW), colorrange = (0, 1))
Colorbar(f[1, 2], hm; label = "Model score")
lines!(ax, amerinord; color = :grey10)
hidedecorations!(ax)
hidespines!(ax)
f

# Vulpes vulpes
    #Import GBIF occurences
        records2 = GBIF.download("10.15468/dl.ggsnrw")
        renard = taxon(first(entity(records2)))

            #Get the Phylopic silouhette
            silhouette_renard = Phylopic.imagesof(renard; items = 1)
            Phylopic.attribution(silhouette_renard)

            #Visualisation
            f = Figure(; size = (600, 300))
            ax = Axis(f[1, 1]; aspect = DataAspect())
            poly!(ax, amerinord; color = :grey95)
            scatter!(ax, records2; color = :forestgreen, label = "Occurences V. vulpes")
            hidedecorations!(ax)
            hidespines!(ax)
            lines!(ax, amerinord; color = :grey10)
            axislegend(ax; position = :rt, framevisible = false)
            silhouetteplot!(
                ax,
                extent.left + 0.3,
                extent.top - 0.2,
                silhouette_renard;
                markersize = 70,
                color = :black,
                )
            f

