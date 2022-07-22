# Static
struct StaticVariable{F} <: JuDOVariable{F}
    ℓ::F
    u::F
end

# Dynamic
struct DynamicVariable{F} <: JuDOVariable{F}
    ℓ::F
    u::F
    initial::Tuple{F,F}
    terminal::Tuple{F,F}
end

# Macros