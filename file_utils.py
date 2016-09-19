def lines_to_array(lines):
    return map(lambda x: map(int, x.split(' ')), lines)

def codes_from_file(file_path):
    with open(file_path, 'r') as fp:
        data = fp.read()
        lines = data.split('\n')
 
        start_read_index = 0
        result = dict()
 
        # 4 here because last entry must have at least 5 lines (1 for meta, 1 for
        # s_y_joins, 1 for y_s_joins and 2 blank lines
        while start_read_index < len(lines) - 4:
            R, frame_len, syndrome_len = lines[start_read_index].split(' ')
            R = float(R)
            frame_len = int(frame_len)
            syndrome_len = int(syndrome_len)
 
            begin_s_y_joins_index = start_read_index + 1
            end_s_y_joins_index = begin_s_y_joins_index + syndrome_len
            s_y_joins = lines_to_array(lines[begin_s_y_joins_index:end_s_y_joins_index])
 
            begin_y_s_joins_index = end_s_y_joins_index
            end_y_s_joins_index = begin_y_s_joins_index + frame_len
            y_s_joins = lines_to_array(lines[begin_y_s_joins_index:end_y_s_joins_index])
 
            begin_punct_list_index = end_y_s_joins_index
            end_puct_list_index = begin_punct_list_index + 1
            punct_list = lines_to_array(lines[begin_punct_list_index:end_puct_list_index])
 
            result[(R, frame_len)] = {
                'R': R,
                'frame_len': frame_len,
                'syndrome_len': syndrome_len,
                's_y_joins': s_y_joins,
                'y_s_joins': y_s_joins,
                'punct_list': punct_list[0]
            }
 
            start_read_index = end_puct_list_index + 1
 
        return result	