class Configuration:
    """
    Utility to parse runko configuration files.
    Used in runko similarly to named function arguments.

    Currently supports .ini format from which sections
    io, simulation, grid, problem and particles are read.
    Found items are usable as attributes of Configuration instance.
    """

    _section_names = "io", "simulation", "grid", "problem", "particles"

    def __init__(self, config_path: str):

        self._config_path = config_path

        from configparser import ConfigParser
        import ast

        parser = ConfigParser()
        parser.optionxform = str  # make option names case sensitive

        if not parser.read(config_path):
            raise ValueError(f"No config found: {config_path}")

        for section in self._section_names:
            attributes = {key: ast.literal_eval(value) for key, value in parser.items(section)}
            self.__dict__.update(attributes)

    def __getattr__(self, name):
        """
        This is called when accessed attribute is missing.
        Instead of raising AttributeError this method raises KeyError,
        because it is more descriptive what actually went wrong.

        Attributes starting with __ are treated as proper attributes
        for which AttributeError is raised instead.
        """

        if name.startswith("__"):
            raise AttributeError(f"Configuration does not define: {name}")
        raise KeyError(f"Key '{name}' not found in config: {self._config_path}")
