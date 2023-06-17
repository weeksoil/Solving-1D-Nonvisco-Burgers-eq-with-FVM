% 初期化
clearvars;

% パラメーター
i_max = 100;    % 格子セル数
XL = -1.0;      % 計算領域左端の座標
XR = 1.0;       % 計算領域右端の座標
a = 1.0;        % 線形移流方程式の移流速度
tstop = 0.5;    % 計算停止時刻


% 配列定義
i = 0;                      % セル番号； (i番目セルの左境界の番号はi、右辺境界の番号はi+1とする)
x = zeros(i_max + 1, 1);    % セル境界の座標
u = zeros(i_max, 1);        % セル平均値 （数値解）
ue = zeros(i_max, 1);       % セル平均値 （厳密解）
ul = zeros(i_max + 1, 1);   % セル境界左側の変数値
ur = zeros(i_max + 1, 1);   % セル境界右側の変数値
f = zeros(i_max + 1, 1);    % セル境界の流束
n = 0;                      % 時間ステップ
t = 0;                      % 計算時間


% 線形移流方程式の初期値を選ぶ
sw1 = 1; %"滑らかな分布（正弦波）：１、"不連続な分布（矩形波）：２

dx = (XR - XL) / (i_max - 4.0);% 格子間隔 計算領域外に二つずつ余分なセルを準備。周期的境界条件向けの設定
dt = 0.2 * dx; % 時間刻み

x(1) = XL - 2.0 * dx;% 計算領域外の２セルも考慮した座標を振る。

[x, u] = initc(i_max, x, dx, sw1, u);% 計算格子，時間刻み，初期条件を設定する

ue = exact(sw1, i_max, ue, x, t, dx);%　厳密解を求める

% メインループ
while t <= tstop

    n = n + 1;
    t = t + dt;

    [ul, ur] = reconstruction_pc(i_max, u, ul, ur);% 空間再構築

    f = riemann_roe(i_max, f, ul, ur);  % リーマンソルバー

    u = update(i_max, u, dt, dx, f); % 時間積分

    u = bc(i_max, u); % 境界条件

    ue = exact(sw1, i_max, ue, x, t, dx);%　厳密解を求める

    fprintf("n=%d, t=%f \n", n, t);

    % 動画保存
    if n == 1
        plotconfig(x(1 : end - 1), ue, u,  t)
        v = VideoWriter('u.mp4','MPEG-4');
        v.FrameRate = 40;
        open(v);
    else
        plotconfig(x(1 : end - 1), ue, u,  t)
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

end

% 動画ファイルを閉じる
close(v);

%% 以下ローカル関数

function[f] = f_flux(x)

f = 0.5*x*x;

end


function [x, u] = initc(i_max, x, dx, sw1, u)

for  i = 2 : i_max + 1
    x(i) = x(i - 1) + dx;% 格子点の座標
end

switch  sw1

    case 1  % 初期の変数値（滑らかな分布）

        for  i = 1 : i_max
            u(i) = 0.5*(1.1 + sin(2.* pi *(x(i)-x(3))));% 三番目の要素（XL、XRの座標値）を基準に考える。
        end

    case 2 % 初期の変数値（不連続な分布）

        for  i = 1 : i_max
            u(i) = 0.1;
        end
        for i = i_max / 2 - 10 : i_max / 2 + 10
            u(i) = 1.0;
        end

    otherwise

        disp("初期値を正しく選択されていません。");

end

end

function [ue] = exact(sw1, i_max, ue, x, t, dx)

switch sw1

    case 1 % 初期の変数値（滑らかな分布）

        for i = 1 : i_max

            c = 2 * pi;
            f = ue(i) - 0.5 * (1.1 + sin(c * (x(i)-x(3) - ue(i) * t)));
            df = 1.0 + 0.5 * c * cos(c *(x(i) - x(3) - ue(i) * t)) * t;
            count = 0;
            while abs(f) >= 1.0e-6
                count = count + 1;
                ue(i) = ue(i) - f / df; % ニュートン法の計算式*/
                f = ue(i) - 0.5 * (1.1 + sin(c*((x(i) - x(3) - ue(i) * t))));
                df = 1.0 + 0.5 * c * cos(c * ((x(i) - x(3) - ue(i) * t))) * t;
                if count > 10000
                    disp("厳密解が収束しないので反復を打ち切ります。");
                    break
                end
            end

        end


    case 2 % 初期の変数値（不連続な分布）

        xc = - dx * 10 + t;
        xl = - dx * 10 + 0.1 * t;
        xr = dx * 10 + 0.55 * t;

        % 周期境界条件　
        if  xl > 1.0
            xl = -2.0 + xl;
        end

        % 周期境界条件　
        if xr > 1.0
            xr = -2.0 + xr;
        end

        % 周期境界条件　
        if xc > 1.0
            xc = -2.0 + xc;
        end

        for i = 1 : i_max
            if x(i) <= xl
                ue(i) = 0.1;
            end
            if x(i) >= xl && x(i) <= xc
                ue(i) =(x(i) - xl) / (xc - xl) * 0.9 + 0.1;
            end
            if x(i) >= xc && x(i) <= xr
                ue(i) = 1.0;
            end
            if x(i) >= xr
                ue(i) = 0.1;
            end
        end

end

end

function [ul, ur] = reconstruction_pc(i_max, u, ul, ur)

for i = 2 : i_max - 2
    ul(i + 1) = u(i); % セル境界(i+1/2)左側の値
    ur(i + 1) = u(i + 1); % セル境界(i+1/2)右側の値
end

end

function [f] = riemann_roe(i_max, f, ul, ur) % 流束を計算する

for i = 3 : i_max - 1
    alpha_12 = 0.5 * (ur(i) + ul(i));
    f(i) = 1.0 / 2.0 * (f_flux(ul(i)) + f_flux(ur(i))) - 1.0 / 2.0 * abs(alpha_12) * (ur(i) - ul(i));
end

end

function [u] = update(i_max, u, dt, dx, f)

for  i = 3 : i_max - 2
    u(i) = u(i) - dt / dx * (f(i + 1) - f(i));% 計算変数を更新する
end

end

function [u] = bc(i_max, u) % 周期境界条件

u(1) = u(i_max - 3); % 計算領域左端の境界条件
u(2) = u(i_max - 2); % 計算領域左端の境界条件
u(i_max - 1) = u(3); % 計算領域右端の境界条件
u(i_max) = u(4); % 計算領域右端の境界条件

end

function [] = plotconfig(x, ue, u, t)

plot(x, ue, x, u)
%plot(x, u)

title(['time = ', num2str(t, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
axis equal;
axis tight;
axis on;
fig=gcf;
fig.Color='white';
ylim([0 1.1]);
xlabel('x')
ylabel('u')

% 凡例
legend({'exact', '1st order Roe'},'Location','southwest','FontSize', 10)

% 新しいプロットの時、軸設定を保持したまま前のグラフィックスオブジェクトを消去
set(gca,'nextplot','replacechildren');

end
