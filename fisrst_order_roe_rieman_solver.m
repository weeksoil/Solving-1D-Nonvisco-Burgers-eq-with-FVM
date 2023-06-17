% ������
clearvars;

% �p�����[�^�[
i_max = 100;    % �i�q�Z����
XL = -1.0;      % �v�Z�̈捶�[�̍��W
XR = 1.0;       % �v�Z�̈�E�[�̍��W
a = 1.0;        % ���`�ڗ��������̈ڗ����x
tstop = 0.5;    % �v�Z��~����


% �z���`
i = 0;                      % �Z���ԍ��G (i�ԖڃZ���̍����E�̔ԍ���i�A�E�Ӌ��E�̔ԍ���i+1�Ƃ���)
x = zeros(i_max + 1, 1);    % �Z�����E�̍��W
u = zeros(i_max, 1);        % �Z�����ϒl �i���l���j
ue = zeros(i_max, 1);       % �Z�����ϒl �i�������j
ul = zeros(i_max + 1, 1);   % �Z�����E�����̕ϐ��l
ur = zeros(i_max + 1, 1);   % �Z�����E�E���̕ϐ��l
f = zeros(i_max + 1, 1);    % �Z�����E�̗���
n = 0;                      % ���ԃX�e�b�v
t = 0;                      % �v�Z����


% ���`�ڗ��������̏����l��I��
sw1 = 1; %"���炩�ȕ��z�i�����g�j�F�P�A"�s�A���ȕ��z�i��`�g�j�F�Q

dx = (XR - XL) / (i_max - 4.0);% �i�q�Ԋu �v�Z�̈�O�ɓ���]���ȃZ���������B�����I���E���������̐ݒ�
dt = 0.2 * dx; % ���ԍ���

x(1) = XL - 2.0 * dx;% �v�Z�̈�O�̂Q�Z�����l���������W��U��B

[x, u] = initc(i_max, x, dx, sw1, u);% �v�Z�i�q�C���ԍ��݁C����������ݒ肷��

ue = exact(sw1, i_max, ue, x, t, dx);%�@�����������߂�

% ���C�����[�v
while t <= tstop

    n = n + 1;
    t = t + dt;

    [ul, ur] = reconstruction_pc(i_max, u, ul, ur);% ��ԍč\�z

    f = riemann_roe(i_max, f, ul, ur);  % ���[�}���\���o�[

    u = update(i_max, u, dt, dx, f); % ���Ԑϕ�

    u = bc(i_max, u); % ���E����

    ue = exact(sw1, i_max, ue, x, t, dx);%�@�����������߂�

    fprintf("n=%d, t=%f \n", n, t);

    % ����ۑ�
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

% ����t�@�C�������
close(v);

%% �ȉ����[�J���֐�

function[f] = f_flux(x)

f = 0.5*x*x;

end


function [x, u] = initc(i_max, x, dx, sw1, u)

for  i = 2 : i_max + 1
    x(i) = x(i - 1) + dx;% �i�q�_�̍��W
end

switch  sw1

    case 1  % �����̕ϐ��l�i���炩�ȕ��z�j

        for  i = 1 : i_max
            u(i) = 0.5*(1.1 + sin(2.* pi *(x(i)-x(3))));% �O�Ԗڂ̗v�f�iXL�AXR�̍��W�l�j����ɍl����B
        end

    case 2 % �����̕ϐ��l�i�s�A���ȕ��z�j

        for  i = 1 : i_max
            u(i) = 0.1;
        end
        for i = i_max / 2 - 10 : i_max / 2 + 10
            u(i) = 1.0;
        end

    otherwise

        disp("�����l�𐳂����I������Ă��܂���B");

end

end

function [ue] = exact(sw1, i_max, ue, x, t, dx)

switch sw1

    case 1 % �����̕ϐ��l�i���炩�ȕ��z�j

        for i = 1 : i_max

            c = 2 * pi;
            f = ue(i) - 0.5 * (1.1 + sin(c * (x(i)-x(3) - ue(i) * t)));
            df = 1.0 + 0.5 * c * cos(c *(x(i) - x(3) - ue(i) * t)) * t;
            count = 0;
            while abs(f) >= 1.0e-6
                count = count + 1;
                ue(i) = ue(i) - f / df; % �j���[�g���@�̌v�Z��*/
                f = ue(i) - 0.5 * (1.1 + sin(c*((x(i) - x(3) - ue(i) * t))));
                df = 1.0 + 0.5 * c * cos(c * ((x(i) - x(3) - ue(i) * t))) * t;
                if count > 10000
                    disp("���������������Ȃ��̂Ŕ�����ł��؂�܂��B");
                    break
                end
            end

        end


    case 2 % �����̕ϐ��l�i�s�A���ȕ��z�j

        xc = - dx * 10 + t;
        xl = - dx * 10 + 0.1 * t;
        xr = dx * 10 + 0.55 * t;

        % �������E�����@
        if  xl > 1.0
            xl = -2.0 + xl;
        end

        % �������E�����@
        if xr > 1.0
            xr = -2.0 + xr;
        end

        % �������E�����@
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
    ul(i + 1) = u(i); % �Z�����E(i+1/2)�����̒l
    ur(i + 1) = u(i + 1); % �Z�����E(i+1/2)�E���̒l
end

end

function [f] = riemann_roe(i_max, f, ul, ur) % �������v�Z����

for i = 3 : i_max - 1
    alpha_12 = 0.5 * (ur(i) + ul(i));
    f(i) = 1.0 / 2.0 * (f_flux(ul(i)) + f_flux(ur(i))) - 1.0 / 2.0 * abs(alpha_12) * (ur(i) - ul(i));
end

end

function [u] = update(i_max, u, dt, dx, f)

for  i = 3 : i_max - 2
    u(i) = u(i) - dt / dx * (f(i + 1) - f(i));% �v�Z�ϐ����X�V����
end

end

function [u] = bc(i_max, u) % �������E����

u(1) = u(i_max - 3); % �v�Z�̈捶�[�̋��E����
u(2) = u(i_max - 2); % �v�Z�̈捶�[�̋��E����
u(i_max - 1) = u(3); % �v�Z�̈�E�[�̋��E����
u(i_max) = u(4); % �v�Z�̈�E�[�̋��E����

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

% �}��
legend({'exact', '1st order Roe'},'Location','southwest','FontSize', 10)

% �V�����v���b�g�̎��A���ݒ��ێ������܂ܑO�̃O���t�B�b�N�X�I�u�W�F�N�g������
set(gca,'nextplot','replacechildren');

end
