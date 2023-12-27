from pathlib import Path
import logging
import logging.config
import yaml

def load_config(path: str) -> dict:
    path_config = Path(path)
    with open(path_config, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.FullLoader)
    return config


def setup_logger(name: str, level="DEBUG", config_path=None) -> logging.Logger:
    """
    Create a logger with a given name
    Args:
        name: (str) logger name
        level: (str) logger level (default DEBUG)
        config_path: (str) Default None. YAML logger configuration file. If no configuration is provided,
                    a basic default is used.
    Returns:
        logger
    """

    logger = logging.getLogger(name)

    if config_path is not None:
        logging_config = load_config(config_path)
        logging.config.dictConfig(logging_config)
    else:
        kwargs = {
            "format": "%(asctime)s - %(name)-7s - %(levelname)s - %(message)s",
            "datefmt": "%Y-%m-%d %H:%M:%S"
        }

        logging.basicConfig(**kwargs)

    logger.setLevel(level)
    return logger


def get_query_folder() -> Path:
    cf = Path(__file__).parent.absolute()
    qdir = cf / 'queries'
    return qdir


def get_query_from_sql_file(
        file_name: str,
        kwargs: dict = None
) -> str:
    with open(get_query_folder() / file_name, 'r') as f:
        query = f.read()
    # injecting kwargs if not None
    if kwargs is not None:
        query = query.format(**kwargs)
    return query
