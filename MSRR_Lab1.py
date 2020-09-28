import plotly.graph_objects as go
from plotly.subplots import make_subplots
from math import sqrt, log2, log10
from scipy.stats import uniform
import numpy as np

def algorithms_fig(means):
    """
    means = dict(sum_eb=[],
                 sum_mt=[],
                 sum_pf=[],
                 min_eb=[],
                 min_mt=[],
                 min_pf=[],
                 avg_eb=[],
                 avg_mt=[],
                 avg_pf=[])
    :param means:
    :return: None
    """
    x = [2 ** i for i in range(7)]

    fig_sum = go.Figure()
    fig_sum.add_trace(go.Scatter(x=x, y=means['sum_eb'], name='Sum EB'))
    fig_sum.add_trace(go.Scatter(x=x, y=means['sum_mt'], name='Sum MT'))
    fig_sum.add_trace(go.Scatter(x=x, y=means['sum_pf'], name='Sum PF'))
    fig_sum.update_layout(title_text='D_sum', width=1000, height=800)
    fig_sum.update_xaxes(title_text='N', type='log')
    fig_sum.update_yaxes(title_text='D_sum, Мбит/с')
    fig_sum.show()

    fig_min = go.Figure()
    fig_min.add_trace(go.Scatter(x=x, y=means['min_eb'], name='Min EB'))
    fig_min.add_trace(go.Scatter(x=x, y=means['min_mt'], name='Min MT'))
    fig_min.add_trace(go.Scatter(x=x, y=means['min_pf'], name='Min PF'))
    fig_min.update_layout(title_text='D_min', width=1000, height=800)
    fig_min.update_xaxes(title_text='N', type='log')
    fig_min.update_yaxes(title_text='D_min, Мбит/с')
    fig_min.show()

    fig_avg = go.Figure()
    fig_avg.add_trace(go.Scatter(x=x, y=means['avg_eb'], name='Avg EB'))
    fig_avg.add_trace(go.Scatter(x=x, y=means['avg_mt'], name='Avg MT'))
    fig_avg.add_trace(go.Scatter(x=x, y=means['avg_pf'], name='Avg PF'))
    fig_avg.update_layout(title_text='D_avg', width=1000, height=800)
    fig_avg.update_xaxes(title_text='N', type='log')
    fig_avg.update_yaxes(title_text='D_avg, Мбит/с')
    fig_avg.show()

def distribution_features_fig(sample: tuple, max_radius: int):
    """
    Гистограмма и эмпирическая функция распределения
    :param sample: Выборка координат
    :param max_radius: Максимальный радиус
    """

    x = [i[1] for i in sample]
    true_x = [i for i in range(2001)]
    true_pdf = [2 * i / max_radius ** 2 for i in true_x]
    true_cdf = [i ** 2 / max_radius ** 2 for i in true_x]

    fig = make_subplots(rows=1,
                        cols=2,
                        subplot_titles=("PDF",
                                        "CDF"))

    hist_trace = go.Histogram(x=x,
                              histnorm='probability density',
                              name='Probability Density Histogram',
                              xbins=dict(start=0,
                                         end=max_radius,
                                         size=max_radius / 10))

    true_pdf_trace = go.Scatter(x=true_x,
                                y=true_pdf,
                                name='True PDF')

    empirical_df_trace = go.Histogram(x=x,
                                      histnorm='probability density',
                                      name='Empirical CDF',
                                      cumulative_enabled=True)

    cdf_trace = go.Scatter(x=true_x,
                           y=true_cdf,
                           name='True CDF')

    fig.append_trace(hist_trace, 1, 1)
    fig.append_trace(true_pdf_trace, 1, 1)
    fig.append_trace(empirical_df_trace, 1, 2)
    fig.append_trace(cdf_trace, 1, 2)

    fig.update_layout(title_text='Distribution Features')
    fig.show()

def a(constants, test):
    """
    result = (1.1 * lg(f_0) - 0.7) * h_rx - (1.56 * lg(f_0) - 0.8
    :return: result
    """
    result = (1.1 * log10(constants['f_0']) - 0.7) * \
             constants['h_rx'] - \
             (1.56 * log10(constants['f_0']) - 0.8)

    if test:
        print("a = {}".format(result))

    return result

def lg_lx10(constants, radius, test):
    """
    result = 46.3 + 33.9 * lf(f_0) - 13.82 * lg(h_bs) - a(h_rx) +
    + (44.9 - 6.55 * lg(h_rx)) * lg(d) + s
    s = 3
    d = sqrt(radius ** 2 + (h_bs - h_rx) ** 2)
    :return: result
    """
    result = 46.3 + 33.9 * log10(constants['f_0']) - 13.82 * \
             log10(constants['h_bs']) - a(constants, test) + \
             (44.9 - 6.55 * log10(constants['h_rx'])) * \
             log10(sqrt((radius ** 2 +
                         (constants['h_bs'] - constants['h_rx']) ** 2)
                        / 1000)) + constants['s']

    if test:
        print("10 * lg(L) = {}\n".format(result))

    return result

def p_rx(constants, radius, test):
    """
    result = P_TX / L
    P_TX = 320
    L = 10 ** (10 * lg(L) / 10)
    """

    result = constants['p_tx'] / \
             (10 ** (lg_lx10(constants, radius, test) / 10))

    if test:
        print("P_RX = {}".format(result))

    return result

def p_n(constants, test):
    """
    result = df * T * k * k_n
    df =
    T = 300
    k = 1.38e-23
    k_n = 3
    :return: result
    """

    result = constants['df'] * constants['t'] * \
             constants['k'] * constants['k_n']

    if test:
        print("P_N = {}".format(result))

    return result

def snr(constants, radius, test):
    """
    result = p_rx / p_n
    :return: result
    """

    result = p_rx(constants, radius, test) / p_n(constants, test)

    if test:
        print("SNR = {}\n".format(result))

    return result

def channel_capacity(constants, radius, test=False):
    """
    lg_lx10

    C = df * log2(1 + SNR)
    """

    result = constants['df'] * log2(1 + snr(constants, radius, test))

    if test:
        print("C = {0}\n".format(result))

    return result

# Генерация координат клиентов
def coord_generator(low, high, count, test=False):
    """
    :param low: нижняя граница генерации
    :param high: верхняя граница геенрации
    :param count: количество координат
    :param test: вывод примера
    :return: tuple(array([x, y]), radius)
    """

    generated = 0
    coords = []

    while generated < count:
        x_y = uniform.rvs(loc=low, scale=high-low, size=2)  # Координаты
        radius = sqrt(x_y[0] ** 2 + x_y[1] ** 2)
        if radius <= 2000:
            generated += 1
            x_y_r = tuple([x_y, radius])
            coords += [x_y_r]

    if test:
        print("Координаты: \nx = {0}, y = {1}, R = {2}".format(
            coords[0][0][0], coords[0][0][1], coords[0][1]
        ))

    return tuple(coords)

# Пример распределения АБ в радиусе
def placement_example_fig(constants):

    radius = constants['radius']
    circle_x = [i - 2000 for i in range(radius - (-radius))]
    rev_circle_x = circle_x.copy()
    rev_circle_x.reverse()
    upper_half_circle_y = [sqrt(radius ** 2 - x ** 2)
                           for x in circle_x]
    lower_half_circle_y = [-sqrt(radius ** 2 - x ** 2)
                           for x in rev_circle_x]
    circle_x = circle_x + rev_circle_x
    circle_y = upper_half_circle_y + lower_half_circle_y

    clients = coord_generator(-radius, radius, 2 ** 6)

    clients_x = [client[0][0] for client in clients]
    clients_y = [client[0][1] for client in clients]

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=circle_x, y=circle_y,
                             name="Радиус вокруг вышки", mode="lines"))
    fig.add_trace(go.Scatter(x=clients_x, y=clients_y,
                             name="Абоненты", mode="markers"))
    fig.add_trace(go.Scatter(x=[0], y=[0],
                             name="Базовая станция", mode="markers"))

    fig.update_layout(title_text="Пример распределения абонентов",
                      width=1000, height=800)
    fig.update_xaxes(range=[-2 * radius, 2 * radius])
    fig.update_yaxes(range=[-2 * radius, 2 * radius])

    fig.show()

# Тест генерации на соответствие нужному распределению
def test_distribution(radius):
    test_coords = coord_generator(-radius, radius, 2 ** 11)
    distribution_features_fig(test_coords, radius)

def equal_blind(c_clients):
    d = [1 / sum([1 / c for c in c_clients]) for i in c_clients]
    return d

def maximum_throughput(c_clients):
    d = [max(c_clients)] + [0 for i in range(len(c_clients) - 1)]
    return d

def proportion_fair(c_clients):
    d = [c / len(c_clients) for c in c_clients]
    return d

# Моделирование
def algorithms_calculation(constants):
    means = dict(sum_eb=[],
                 sum_mt=[],
                 sum_pf=[],
                 min_eb=[],
                 min_mt=[],
                 min_pf=[],
                 avg_eb=[],
                 avg_mt=[],
                 avg_pf=[])
    for i in range(7):
        d = dict(sum_eb=[],
                 sum_mt=[],
                 sum_pf=[],
                 min_eb=[],
                 min_mt=[],
                 min_pf=[],
                 avg_eb=[],
                 avg_mt=[],
                 avg_pf=[])
        for j in range(2 ** 10):
            clients = coord_generator(-constants['radius'],
                                      constants['radius'], 2 ** i)
            c_clients = [channel_capacity(constants, client[1])
                         for client in clients]
            d_eb = equal_blind(c_clients)
            d_mt = maximum_throughput(c_clients)
            d_pf = proportion_fair(c_clients)
            d['sum_eb'] += [np.sum(d_eb)]
            d['sum_mt'] += [np.sum(d_mt)]
            d['sum_pf'] += [np.sum(d_pf)]
            d['min_eb'] += [np.min(d_eb)]
            d['min_mt'] += [np.min(d_mt)]
            d['min_pf'] += [np.min(d_pf)]
            d['avg_eb'] += [np.mean(d_eb)]
            d['avg_mt'] += [np.mean(d_mt)]
            d['avg_pf'] += [np.mean(d_pf)]
        for key in means.keys():
            means[key] += [np.mean(d[key])]
    return means

# Пример расчетов
def calculation_test(constants):
    client = coord_generator(
        -constants['radius'], constants['radius'], 1, True
    )
    channel_capacity(constants, client[0][1], True)


if __name__ == "__main__":
    """
    Окумура-Хата (Large City)
    """
    consts = dict(radius=2000,
                     p_tx=320,
                     f_0=1800,
                     df=5,
                     k_n=3,
                     k=1.38e-23,
                     h_bs=30,
                     h_rx=1.5,
                     t=300,
                     s=3)

    test_distribution(consts['radius'])
    calculation_test(consts)

    placement_example_fig(consts)

    results = algorithms_calculation(consts)
    algorithms_fig(results)

