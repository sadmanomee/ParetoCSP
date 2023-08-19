import hyperopt as hy


def hy_parameter_setting(label, config, ptype='float'):
    if type(config) is list:
        if ptype == 'int':
            parameter = hy.hp.choice(label, list(range(config[0], config[1] + 1)))
        else:
            parameter = hy.hp.uniform(label, config[0], config[1])
    elif type(config) is tuple:
        parameter = hy.hp.choice(label, config)
    elif type(config) is int or type(config) is float:
        parameter = hy.hp.choice(label, [config, ])
    else:
        raise Exception("Parameter `" + label + "` setting error! Please check!")
    return parameter


def sko_parameter_setting():
    pass
