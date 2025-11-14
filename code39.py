import sys
import json

sys.setrecursionlimit(1000000)

def check_hamming_distance_one(g1: str, g2: str) -> bool:
    """
    Проверяет, равно ли расстояние Хэмминга между двумя строками ровно 1.
    Предполагается, что строки имеют одинаковую длину.
    """
    diff_count = 0
    for i in range(len(g1)):
        if g1[i] != g2[i]:
            diff_count += 1
            if diff_count > 1:
                # Нашли второе различие, нет смысла проверять дальше
                return False
                
    return diff_count == 1

def build_h1_graph(genotypes: list[str]) -> dict[int, list[int]]:
    """
    Строит граф H=1 (расстояние Хэмминга 1) для списка генотипов.
    
    Возвращает список смежности (dict), где ключ - индекс генотипа,
    а значение - список индексов его соседей (H=1).
    """
    n = len(genotypes)
    adjacency_list = {i: [] for i in range(n)}
    
    if n == 0:
        return {}

    print(f"Начинаем построение графа для {n} узлов...")
    
    for i in range(n):
        # Логирование прогресса (полезно для больших N)
        if n > 1000 and i % 1000 == 0:
            print(f"Обработка узла {i} / {n}")
            
        for j in range(i + 1, n):
            if check_hamming_distance_one(genotypes[i], genotypes[j]):
                adjacency_list[i].append(j)
                adjacency_list[j].append(i)
                
    print("Построение графа завершено.")
    return adjacency_list

# --- НОВАЯ ФУНКЦИЯ ---
def find_connected_components(graph: dict[int, list[int]], n_nodes: int) -> list[list[int]]:
    """
    Находит все компоненты связности в графе (заданном списком смежности).
    Использует Поиск в глубину (DFS).
    
    Возвращает список списков, где каждый внутренний список - 
    это индексы узлов в одной компоненте связности.
    """
    visited = set()
    all_components = []

    # Рекурсивная функция DFS для обхода одной компоненты
    def dfs(node: int, current_component: list[int]):
        visited.add(node)
        current_component.append(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, current_component)

    # Проходим по всем узлам (генотипам) от 0 до N-1
    for i in range(n_nodes):
        if i not in visited:
            # Нашли узел, который еще не посещали -> это начало новой компоненты
            new_component = []
            dfs(i, new_component)
            all_components.append(new_component)
    
    return all_components

def main():
    if len(sys.argv) < 2:
        print("Usage: python find_max_subset.py <landscape_file>")
        sys.exit(1)

    fname = sys.argv[1]
    genotypes = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            genotypes.append(line.split()[0])

    # 1. Строим граф
    h1_graph = build_h1_graph(genotypes)

    # 2. Находим компоненты связности
    components = find_connected_components(h1_graph, len(genotypes))

    print("\n--- Анализ компонент связности ---")

    # 1) Количество компонент связности
    print(f"\n1. Общее количество компонент связности: {len(components)}")

    # 2) Количество генотипов в каждой из компонент связности
    component_sizes = [len(comp) for comp in components]
    print(f"\n2. Количество генотипов в компонентах: {component_sizes}")

    # 3) Список индексов генотипов, входящих в каждую компоненту связности
    print("\n3. Cписки индексов генотипов по компонентам:")
    for i, comp in enumerate(components):
        # Сортируем для более чистого вывода
        comp.sort() 
        print(f"   Компонента {i+1} (размер {len(comp)}): {comp}")

    # Ожидаемый вывод для этого примера:
    #   Компонента 1 (размер 5): [0, 1, 2, 3, 4]
    #   Компонента 2 (размер 1): [5]

if __name__ == "__main__":
    main()
