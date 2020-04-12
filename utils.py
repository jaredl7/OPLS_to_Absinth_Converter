import os
import yaml


def validate_config_file(yaml_config_filename):
    required = 'opls-topology-file,opls-parameters-file,abs-opls-pka-file,abs-charmm-pka-file'.split(',')
    parsed_config = load_yaml_file(yaml_config_filename)
    missing_params = list(set(required).difference(list(parsed_config.keys())))
    if len(missing_params) > 0:
        raise RuntimeError('Missing parameters: %s. Exiting.' % '; '.join(missing_params))
    
    for param in parsed_config:
        param_value = parsed_config[param]
        if param_value is None:
            raise RuntimeError('The value for paramter "{}" is empty. Exiting.'.format(param))
    return parsed_config


def load_yaml_file(filename):
    with open(filename, 'r') as config_file:
        data = yaml.load(config_file)
    return data


def load_config(yaml_config_filename):
    return validate_config_file(yaml_config_filename)


def inline_comment_location(line, inline_comment='!'):
    index = None
    if inline_comment in line:
        index = line.index(inline_comment)
    return index


def strip_comments(filename):
    text = open(filename, 'r').read()
    data = text.split('\n')

    # remove the lines with comments
    data_no_comments = list()
    comments_positions = list()
    for line_index, raw_line in enumerate(data):
        line = raw_line.strip()

        if line.startswith('*') or line.startswith('!') or line.startswith('#'):
            comments_positions.append(line_index)
        else:
            # some of the lines contain in-line comments. Remove them
            line_no_comment = line[::]
            if '!' in line:
                inline_comment_index = inline_comment_location(line, inline_comment='!')
                line_no_comment = line[:inline_comment_index].strip()
            if '#' in line_no_comment:
                inline_comment_index = inline_comment_location(line_no_comment, inline_comment='#')
                line_no_comment = line_no_comment[:inline_comment_index].strip()
            data_no_comments.append(line_no_comment)

    # remove the line containing "END"
    if 'END' in data_no_comments[-1]:
        data_no_comments.pop()
    return data, data_no_comments, comments_positions