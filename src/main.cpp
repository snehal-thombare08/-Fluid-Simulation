#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <cstdlib>

const int WIDTH = 1200, HEIGHT = 800;
const int N = 120;
const float DX = (float)WIDTH / N;
const float DY = (float)HEIGHT / N;
const float DT = 0.1f;
const float VISC = 0.00002f;
const float DIFF = 0.00001f;
const int ITER = 8;

struct Fluid {
    int n;
    float dt, diff, visc;
    std::vector<float> s, dens;
    std::vector<float> vx, vy, vx0, vy0;

    Fluid(int n, float dt, float diff, float visc)
        : n(n), dt(dt), diff(diff), visc(visc),
          s(n*n,0), dens(n*n,0),
          vx(n*n,0), vy(n*n,0), vx0(n*n,0), vy0(n*n,0) {}

    int IX(int x, int y) { return x + y*n; }

    void setBnd(int b, std::vector<float>& x) {
        for(int i=1;i<n-1;i++){
            x[IX(i,0)]   = b==2 ? -x[IX(i,1)]   : x[IX(i,1)];
            x[IX(i,n-1)] = b==2 ? -x[IX(i,n-2)] : x[IX(i,n-2)];
            x[IX(0,i)]   = b==1 ? -x[IX(1,i)]   : x[IX(1,i)];
            x[IX(n-1,i)] = b==1 ? -x[IX(n-2,i)] : x[IX(n-2,i)];
        }
        x[IX(0,0)]     = 0.5f*(x[IX(1,0)]+x[IX(0,1)]);
        x[IX(0,n-1)]   = 0.5f*(x[IX(1,n-1)]+x[IX(0,n-2)]);
        x[IX(n-1,0)]   = 0.5f*(x[IX(n-2,0)]+x[IX(n-1,1)]);
        x[IX(n-1,n-1)] = 0.5f*(x[IX(n-2,n-1)]+x[IX(n-1,n-2)]);
    }

    void linSolve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c) {
        float cRecip = 1.f/c;
        for(int k=0;k<ITER;k++){
            for(int j=1;j<n-1;j++)
                for(int i=1;i<n-1;i++)
                    x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i+1,j)]+x[IX(i-1,j)]+x[IX(i,j+1)]+x[IX(i,j-1)]))*cRecip;
            setBnd(b,x);
        }
    }

    void diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff2) {
        float a = dt*diff2*(n-2)*(n-2);
        linSolve(b,x,x0,a,1+6*a);
    }

    void advect(int b, std::vector<float>& d, std::vector<float>& d0,
                std::vector<float>& vx2, std::vector<float>& vy2) {
        float dtx=dt*(n-2), dty=dt*(n-2);
        for(int j=1;j<n-1;j++) for(int i=1;i<n-1;i++){
            float x=i-dtx*vx2[IX(i,j)], y=j-dty*vy2[IX(i,j)];
            x=std::max(0.5f,std::min((float)(n-1)-0.5f,x));
            y=std::max(0.5f,std::min((float)(n-1)-0.5f,y));
            int i0=(int)x, i1=i0+1, j0=(int)y, j1=j0+1;
            float sx=x-i0, sy=y-j0;
            d[IX(i,j)] = (1-sx)*((1-sy)*d0[IX(i0,j0)]+sy*d0[IX(i0,j1)])
                        +    sx *((1-sy)*d0[IX(i1,j0)]+sy*d0[IX(i1,j1)]);
        }
        setBnd(b,d);
    }

    void project(std::vector<float>& vx2, std::vector<float>& vy2,
                 std::vector<float>& p,   std::vector<float>& div) {
        float h = 1.f/(n-2);
        for(int j=1;j<n-1;j++) for(int i=1;i<n-1;i++){
            div[IX(i,j)] = -0.5f*h*(vx2[IX(i+1,j)]-vx2[IX(i-1,j)]+vy2[IX(i,j+1)]-vy2[IX(i,j-1)]);
            p[IX(i,j)]=0;
        }
        setBnd(0,div); setBnd(0,p);
        linSolve(0,p,div,1,6);
        for(int j=1;j<n-1;j++) for(int i=1;i<n-1;i++){
            vx2[IX(i,j)] -= 0.5f*(p[IX(i+1,j)]-p[IX(i-1,j)])/(h*(n-2));
            vy2[IX(i,j)] -= 0.5f*(p[IX(i,j+1)]-p[IX(i,j-1)])/(h*(n-2));
        }
        setBnd(1,vx2); setBnd(2,vy2);
    }

    void step() {
        diffuse(1,vx0,vx,visc); diffuse(2,vy0,vy,visc);
        project(vx0,vy0,vx,vy);
        advect(1,vx,vx0,vx0,vy0); advect(2,vy,vy0,vx0,vy0);
        project(vx,vy,vx0,vy0);
        diffuse(0,s,dens,diff); advect(0,dens,s,vx,vy);
    }

    void addDensity(int x, int y, float amount) {
        if(x<0||x>=n||y<0||y>=n) return;
        dens[IX(x,y)] += amount;
    }

    void addVelocity(int x, int y, float ax, float ay) {
        if(x<0||x>=n||y<0||y>=n) return;
        vx[IX(x,y)] += ax; vy[IX(x,y)] += ay;
    }

    void fadeDensity() {
        for(auto& d : dens) d = std::max(0.f, d*0.995f);
    }
};

sf::Color fluidColor(float d, int mode) {
    float t = std::min(1.f, d * 0.15f);
    if(mode==0) {
        return sf::Color((uint8_t)(t*80),(uint8_t)(t*180),(uint8_t)(100+t*155));
    } else if(mode==1) {
        return sf::Color((uint8_t)(t*255),(uint8_t)(t*80+20),(uint8_t)(t*20));
    } else if(mode==2) {
        return sf::Color((uint8_t)(t*255),(uint8_t)(t*255),(uint8_t)(t*255));
    } else {
        return sf::Color(
            (uint8_t)(127+127*std::sin(t*6.28f)),
            (uint8_t)(127+127*std::sin(t*6.28f+2.09f)),
            (uint8_t)(127+127*std::sin(t*6.28f+4.19f))
        );
    }
}

int main(){
    sf::RenderWindow window(
        sf::VideoMode({(unsigned)WIDTH,(unsigned)HEIGHT}),
        "Fluid Simulation | LMB:Add Fluid | RMB:Erase | 1-4:Color | V:Velocity | Space:Clear"
    );
    window.setFramerateLimit(60);

    Fluid fluid(N, DT, DIFF, VISC);

    std::vector<uint8_t> pixels(WIDTH*HEIGHT*4, 0);
    sf::Texture tex({(unsigned)WIDTH,(unsigned)HEIGHT});
    sf::Sprite sprite(tex);

    int colorMode = 0;
    bool showVel = false;
    sf::Vector2i prevMp{0,0};

    sf::Font font; bool hf = font.openFromFile("C:/Windows/Fonts/arial.ttf");
    sf::Text hud(font); hud.setCharacterSize(14);
    hud.setFillColor(sf::Color::White);
    hud.setOutlineColor(sf::Color::Black); hud.setOutlineThickness(1.5f);
    hud.setPosition({8.f,8.f});

    sf::Clock clock;
    bool lmb=false, rmb=false;

    while(window.isOpen()){
        clock.restart();
        auto mp = sf::Mouse::getPosition(window);
        int fx = (int)(mp.x/DX), fy = (int)(mp.y/DY);

        while(auto ev = window.pollEvent()){
            if(ev->is<sf::Event::Closed>()) window.close();
            if(auto* k = ev->getIf<sf::Event::KeyPressed>()){
                if(k->code==sf::Keyboard::Key::Escape) window.close();
                if(k->code==sf::Keyboard::Key::Space){
                    std::fill(fluid.dens.begin(),fluid.dens.end(),0.f);
                    std::fill(fluid.vx.begin(),fluid.vx.end(),0.f);
                    std::fill(fluid.vy.begin(),fluid.vy.end(),0.f);
                }
                if(k->code==sf::Keyboard::Key::Num1) colorMode=0;
                if(k->code==sf::Keyboard::Key::Num2) colorMode=1;
                if(k->code==sf::Keyboard::Key::Num3) colorMode=2;
                if(k->code==sf::Keyboard::Key::Num4) colorMode=3;
                if(k->code==sf::Keyboard::Key::V) showVel=!showVel;
            }
            if(auto* mb=ev->getIf<sf::Event::MouseButtonPressed>()){
                if(mb->button==sf::Mouse::Button::Left) lmb=true;
                if(mb->button==sf::Mouse::Button::Right) rmb=true;
            }
            if(auto* mb=ev->getIf<sf::Event::MouseButtonReleased>()){
                if(mb->button==sf::Mouse::Button::Left) lmb=false;
                if(mb->button==sf::Mouse::Button::Right) rmb=false;
            }
        }

        if(lmb && fx>0 && fx<N-1 && fy>0 && fy<N-1){
            float dvx = (float)(mp.x - prevMp.x)*1.5f;
            float dvy = (float)(mp.y - prevMp.y)*1.5f;
            for(int dy=-2;dy<=2;dy++) for(int dx=-2;dx<=2;dx++){
                fluid.addDensity(fx+dx,fy+dy,500.f);
                fluid.addVelocity(fx+dx,fy+dy,dvx,dvy);
            }
        }
        if(rmb && fx>0 && fx<N-1 && fy>0 && fy<N-1){
            for(int dy=-3;dy<=3;dy++) for(int dx=-3;dx<=3;dx++){
                int cx=fx+dx, cy=fy+dy;
                if(cx>=0&&cx<N&&cy>=0&&cy<N)
                    fluid.dens[fluid.IX(cx,cy)]=0;
            }
        }

        prevMp = mp;
        fluid.step();
        fluid.fadeDensity();

        // Render
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            sf::Color col = fluidColor(fluid.dens[fluid.IX(i,j)], colorMode);
            int px=(int)(i*DX), py=(int)(j*DY);
            int pw=(int)DX+1, ph=(int)DY+1;
            for(int dy=0;dy<ph&&py+dy<HEIGHT;dy++)
                for(int dx=0;dx<pw&&px+dx<WIDTH;dx++){
                    int idx=((py+dy)*WIDTH+(px+dx))*4;
                    pixels[idx]=col.r; pixels[idx+1]=col.g;
                    pixels[idx+2]=col.b; pixels[idx+3]=255;
                }
        }
        tex.update(pixels.data());

        window.clear(sf::Color(8,10,18));
        window.draw(sprite);

        // Velocity field overlay
        if(showVel){
            sf::VertexArray lines(sf::PrimitiveType::Lines);
            for(int j=2;j<N-2;j+=3) for(int i=2;i<N-2;i+=3){
                float px2=(i+0.5f)*DX, py2=(j+0.5f)*DY;
                float vx2=fluid.vx[fluid.IX(i,j)]*8.f;
                float vy2=fluid.vy[fluid.IX(i,j)]*8.f;
                lines.append({sf::Vector2f(px2,py2),sf::Color(255,255,255,80)});
                lines.append({sf::Vector2f(px2+vx2,py2+vy2),sf::Color(255,200,50,200)});
            }
            window.draw(lines);
        }

        if(hf){
            hud.setString("LMB: Paint  RMB: Erase  1:Cyan  2:Lava  3:White  4:Rainbow  V:Velocity(" +
                std::string(showVel?"ON":"OFF")+")  Space:Clear");
            window.draw(hud);
        }
        window.display();
    }
    return 0;
}



