import re

def extract_info(path):
    match = re.search(r'(\d{10})/zpc-(\d+)\.root$', path)
    if match:
        date_part = match.group(1)
        number_part = int(match.group(2))
        return date_part, number_part
    return None, None

def main():
    with open('zpc.list', 'r') as file:
        file_paths = file.read().splitlines()
    
    # 提取日期和序号信息，并进行排序
    sorted_paths = sorted(file_paths, key=lambda x: (extract_info(x)[0], extract_info(x)[1]), reverse=True)

    # 将排序后的路径写入 b.txt 文件
    with open('zpc_sorted.list', 'w') as output_file:
        for path in sorted_paths:
            output_file.write(path + '\n')

if __name__ == '__main__':
    main()