#!/bin/bash

# 定义输出文件夹的根目录
base_dir="/storage/fdunphome/wangchunzheng/cz_coal/batch"

# 创建一个目录来存储所有合并后的中间文件
central_dir="${base_dir}/merged_files"
mkdir -p "$central_dir"

# 遍历每个日期文件夹
for date_folder in "$base_dir"/*; do
  if [ -d "$date_folder" ]; then
    # 获取日期文件夹的基本名称，即日期
    date=$(basename "$date_folder")

    # 创建一个针对每个日期的合并命令文件，存放在中心目录中
    merge_script="${central_dir}/merge_commands-${date}.sh"
    : > "$merge_script"  # 清空或创建新文件
    chmod +x "$merge_script"  # 给予执行权限

    # 定义每个日期文件夹中的合并后文件，同样存放在中心目录中
    date_merged_file="${central_dir}/obvs_coalHadrons-${date}.root"

    # 找到所有的观察结果文件
    obs_files=$(find "$date_folder" -name 'obvs_coalHadrons-*.root' -print0 | xargs -0 echo)

    # 如果找到了观察结果文件，生成hadd命令并添加到日期的合并脚本中
    if [ -n "$obs_files" ]; then
      echo "hadd -f $date_merged_file $obs_files" >> "$merge_script"
    fi

    # 运行这个日期的合并脚本
    bash "$merge_script" &
  fi
done

# 等待所有后台合并命令完成
wait

# 现在合并所有的日期文件夹的合并文件
all_merged_files=$(find "$central_dir" -name 'obvs_coalHadrons-*.root')
final_merged_file="${central_dir}/obvs_coalHadrons-FinalMerged.root"

if [ -n "$all_merged_files" ]; then
  hadd -f "$final_merged_file" $all_merged_files
  echo "All observation result files have been merged into $final_merged_file"
else
  echo "No files to merge."
fi