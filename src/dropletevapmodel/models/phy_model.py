from pydantic import BaseModel, Field

class FluidProperties(BaseModel):
    rho_d: float = Field(description="Density", ge=0)
    boiling_t: int = Field(description="Boiling temperature", ge=0)
    cp_l: float = Field(description="Specific heat capacity", ge=0)
    mm_d: int = Field(description="Molecular mass", ge=0)
    sigma_d: float = Field(description="Surface tension", ge=0)
    eps_fact: float = Field(description="Factor", ge=0)
    speed_d: int = Field(description="Speed", ge=0)
    viscosity_a: float = Field(description="Viscosity coefficient a", ge=0)
    viscosity_b: float = Field(description="Viscosity coefficient b", ge=0)
    cp_a: float = Field(description="Heat capacity coefficient a")
    cp_b: float = Field(description="Heat capacity coefficient b")
    cp_c: float = Field(description="Heat capacity coefficient c")
    cp_d: float = Field(description="Heat capacity coefficient d")
    cp_e: float = Field(description="Heat capacity coefficient e")
    k_a: float = Field(description="Thermal conductivity coefficient a")
    k_b: float = Field(description="Thermal conductivity coefficient b")
    k_c: float = Field(description="Thermal conductivity coefficient c")



class GasProperties(BaseModel):
    mm_g: float = Field(description="Molecular mass", ge=0)
    sigma_g: float = Field(description="Surface tension", ge=0)
    eps_fact: float = Field(description="Factor", ge=0)
    r: float = Field(description="Gas constant", ge=0)
    speed_g: int = Field(description="Speed", ge=0)
    viscosity_a: float = Field(description="Viscosity coefficient a", ge=0)
    viscosity_b: float = Field(description="Viscosity coefficient b", ge=0)
    cp_a: float = Field(description="Heat capacity coefficient a")
    cp_b: float = Field(description="Heat capacity coefficient b")
    cp_c: float = Field(description="Heat capacity coefficient c")
    cp_d: float = Field(description="Heat capacity coefficient d")
    cp_e: float = Field(description="Heat capacity coefficient e")
    k_a: float = Field(description="Thermal conductivity coefficient a")
    k_b: float = Field(description="Thermal conductivity coefficient b")
    k_c: float = Field(description="Thermal conductivity coefficient c")

class ModellingConfig(BaseModel):
    timestep: float = Field(description="Simulation timestep", gt=0)
    max_iteration: float = Field(description="Maximum number of iterations", gt=0)

class EnvironmentConfig(BaseModel):
    pressure: float = Field(description="Pressure", gt=0)
    pressure_atm: float = Field(description="Atmospheric pressure", gt=0)
    kb: float = Field(description="Boltzmann constant", gt=0)
    drop_conc: float = Field(description="Droplet concentration", gt=0)

class TZeroConfig(BaseModel):
    d_zero: float = Field(description="Initial value for d_zero", gt=0)
    t_d_zero: float = Field(description="Initial temperature for d_zero", gt=0)
    t_g_zero: float = Field(description="Initial gas temperature", gt=0)
    water_content_g_zero: float = Field(description="Initial water content in gas", ge=0)

class OutputConfig(BaseModel):
    csv: bool = Field(description="Enable CSV output")
    plotting: bool = Field(description="Enable plotting output")


class SimulationConfig(BaseModel):
    modelling: ModellingConfig = Field(description="Modelling configuration")
    environment: EnvironmentConfig = Field(description="Environment configuration")
    t_zero: TZeroConfig = Field(description="Initial conditions")
    output: OutputConfig = Field(description="Output configuration")

class DropletEvapModel(BaseModel):
    fluid_properties: FluidProperties = Field(description="Fluid properties")
    gas_properties: GasProperties = Field(description="Gas properties")
    simulation_config: SimulationConfig = Field(description="Simulation configuration")
