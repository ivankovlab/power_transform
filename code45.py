from collections import deque # Используем deque для эффективной очереди в BFS
import sys

def build_graph(sequences, alphabet):
    """
    Фаза 1: Эффективное построение графа.
    Строит список смежности для графа, где ребра соединяют
    последовательности с расстоянием Хэмминга 1.
    """
    n = len(sequences)
    if n == 0:
        return []
    
    m = len(sequences[0])
    
    # 1. Структуры для быстрого поиска
    # Хеш-таблица (set) для проверки O(1), существует ли строка
    string_set = set(sequences)
    
    # Хеш-карта (dict) для сопоставления строки ее индексу O(1)
    string_to_index = {seq: i for i, seq in enumerate(sequences)}
    
    # 2. Список смежности
    graph = [[] for _ in range(n)]
    
    # Используем set для ребер, чтобы избежать дубликатов 
    # (например, ребра A->B и B->A, добавленного при обходе B)
    edges = set()

    # 3. Построение ребер
    for i, seq in enumerate(sequences):
        original_seq_list = list(seq)
        
        for j in range(m): # Итерация по каждой позиции в строке
            original_char = original_seq_list[j]
            
            for char in alphabet: # Итерация по каждой букве алфавита
                if char == original_char:
                    continue
                
                # Создаем "соседа"
                original_seq_list[j] = char
                neighbor_str = "".join(original_seq_list)
                
                # Проверяем, существует ли этот "сосед" в нашем наборе
                if neighbor_str in string_set:
                    neighbor_index = string_to_index[neighbor_str]
                    
                    # Добавляем ребро. 
                    # Проверяем, что i < neighbor_index, чтобы добавить 
                    # только одно направленное ребро (i, k) и избежать (k, i)
                    if i < neighbor_index:
                        edge = (i, neighbor_index)
                        if edge not in edges:
                            graph[i].append(neighbor_index)
                            graph[neighbor_index].append(i) # Граф неориентированный
                            edges.add(edge)

            # Возвращаем исходный символ на место для следующей итерации
            original_seq_list[j] = original_char
            
    return graph

def find_connected_components(n, graph):
    """
    Фаза 2: Поиск компонент связности с помощью BFS.
    """
    visited = [False] * n
    all_components = []

    for i in range(n):
        if not visited[i]:
            # Нашли узел, который еще не посещали -> начало новой компоненты
            current_component = []
            q = deque([i]) # Очередь для BFS
            visited[i] = True
            
            while q:
                current_node = q.popleft()
                current_component.append(current_node)
                
                # Ищем всех соседей
                for neighbor in graph[current_node]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        q.append(neighbor)
            
            # Сортируем для красивого вывода (опционально)
            current_component.sort()
            all_components.append(current_component)
            
    return all_components

def analyze_sequences(sequences, alphabet):
    """
    Главная функция для анализа и вывода результатов.
    """
    if not sequences:
        print("Список последовательностей пуст.")
        return

    print("Анализ запущен...")
    
    # Фаза 1
    print(f"Построение графа для {len(sequences)} последовательностей...")
    graph = build_graph(sequences, alphabet)
    print("Граф построен.")

    # Фаза 2
    print("Поиск компонент связности...")
    components = find_connected_components(len(sequences), graph)
    print("Поиск завершен.\n")

    # --- Вывод результатов ---
    
    # 1) Число независимых компонент связности
    print(f"## 1. Общее число компонент связности: {len(components)}")
    
    # 2) Количество последовательностей в каждой из компонент
    component_sizes = [len(c) for c in components]
    print(f"\n## 2. Количество последовательностей в компонентах:")
    print(component_sizes)
    
    # 3) Список индексов последовательностей, входящих в каждую компоненту
    print(f"\n## 3. Списки индексов для каждой компоненты:")
    for i, component in enumerate(components):
        print(f"  Компонента {i+1} (размер={len(component)}): {component}")

# Определяем наш алфавит
ALPHABET = "ACDEFGHIKLMNPQRSTVWY*"

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


print("--- ЗАПУСК ---")
analyze_sequences(genotypes, ALPHABET)
