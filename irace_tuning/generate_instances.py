import os


def update_instances_list(instances_dir, instances_list_path):
    """
    Update the instances-list.txt file with all files in the instances directory.

    :param instances_dir: Path to the directory containing instance files.
    :param instances_list_path: Path to the instances-list.txt file.
    """
    # Get all files in the instances directory
    instances = [
        f"instances/{file}"
        for file in os.listdir(instances_dir)
        if file.endswith(".txt")
    ]

    # Write the list of instances to the file
    with open(instances_list_path, "w") as f:
        f.write("# List of instances for irace\n")
        f.write("\n".join(instances))
    print(f"Updated {instances_list_path} with {len(instances)} instances.")


def main():
    base_path = os.getcwd()
    instances_dir = os.path.join(base_path, "instances")
    instances_list_path = os.path.join(base_path, "instances-list.txt")

    update_instances_list(instances_dir, instances_list_path)


if __name__ == "__main__":
    main()
