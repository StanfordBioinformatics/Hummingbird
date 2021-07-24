#!/usr/bin/env python3

from .instance import Instance

MISSING_FIELD_ERR_FMT = "{0} is a required field in the conf file."


class ConfigurationError(Exception):
    def __init__(self, message=None, field=None):
        if not message:
            message = MISSING_FIELD_ERR_FMT.format(field)
        message += ' Please edit your file accordingly and try again.'
        super().__init__(message)


def validate_config_file(conf):
    platform = __verify_present(conf, 'Platform', dict)
    service = __verify_present(platform, 'service', str)

    downsample = conf.get('Downsample')
    if downsample:
        __verify_present(downsample, 'input', dict)
        __verify_present(downsample, 'target', int)
        __verify_present(downsample, 'fractions', list)
        __verify_present(downsample, 'output', str)

    if 'Profiling' in conf:
        profiling_list = __verify_present(conf, 'Profiling', type_check=list)
        for profiling in profiling_list:
            __verify_present(profiling, 'result', str)
            __verify_present(profiling, 'command', str)
            __verify_threads(profiling, service)


def __verify_present(val, key, type_check=None):
    if key not in val:
        raise ConfigurationError(field=key)

    result = val[key]
    if type_check and not isinstance(result, type_check):
        raise ConfigurationError(
            message="{0} field in the conf needs to be of type '{1}'.".format(key, type_check.__name__))
    return result


def __verify_threads(profiling, service):
    threads = __verify_present(profiling, 'thread', list)
    if not Instance.check_threads_supported(threads, service):
        raise ConfigurationError(
            message="Hummingbird does not support threads '{0}' for '{1}' service and only supports '{2}'.".format(
                threads, service, sorted(Instance.get_supported_threads(service))
            ),
        )
