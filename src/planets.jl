# Planet definitions with mean radii in km
# Sources: NASA planetary fact sheets
const Mercury = Planet(:mercury, 2439.7)
const Earth = Planet(:earth, 6371.2)
const Mars = Planet(:mars, 3389.5)
const Jupiter = Planet(:jupiter, 71492.0)
const Saturn = Planet(:saturn, 60268.0)
const Uranus = Planet(:uranus, 25559.0)
const Neptune = Planet(:neptune, 24764.0)
const Ganymede = Planet(:ganymede, 2634.1)

const PLANETS = Dict{Symbol, Planet}(
    :mercury => Mercury,
    :earth => Earth,
    :mars => Mars,
    :jupiter => Jupiter,
    :saturn => Saturn,
    :uranus => Uranus,
    :neptune => Neptune,
    :ganymede => Ganymede,
)

function planet(s::Symbol)
    key = Symbol(lowercase(string(s)))
    haskey(PLANETS, key) || error("Unknown planet: $s. Available: $(keys(PLANETS))")
    return PLANETS[key]
end

planet(p::Planet) = p
