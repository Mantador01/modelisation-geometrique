#include "qte.h"
#include "implicits.h"
#include "ui_interface.h"

#include "blobfield.h"
#include "erosionfield.h"
#include "sphere_tracing.h"

#include "terrain.h"
#include "noise.h"
#include "sdf_nodes.h"

#include <random>         // std::mt19937, std::uniform_real_distribution
#include <QElapsedTimer>  // QElapsedTimer
#include <QDebug>         // qDebug()

#include <QtCore/QElapsedTimer>
#include <QMessageBox>
#include <functional>
#include <climits>

#include <chrono>
#include <iostream>



#include "erosionfield.h"         // batch
#include "erosion_incremental.h"  // incremental

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw->widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->blobsImplicit, SIGNAL(clicked()), this, SLOT(BlobsImplicitExample()));
    connect(uiw->pushButtonErosion, SIGNAL(clicked()), this, SLOT(ErosionImplicitExample()));
    connect(uiw->pushButtonTerrain, SIGNAL(clicked()), this, SLOT(TerrainHeightmapExample()));
    connect(uiw->pushButtonSDF, SIGNAL(clicked()), this, SLOT(SDFExample()));
    connect(uiw->pushButtonBenchmark, SIGNAL(clicked()), this, SLOT(BenchmarkBlobs()));
    connect(uiw->pushButtonErosionIncremental, SIGNAL(clicked()), this, SLOT(RunBenchErosion()));


	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BoxMeshExample()
{
	Mesh boxMesh = Mesh(Box(1.0));

	std::vector<Color> cols;
	cols.resize(boxMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
	UpdateGeometry();
}

void MainWindow::SphereImplicitExample()
{
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::BlobsImplicitExample()
{
//   BlobField implicit;
//   implicit.Clear();
//   implicit.SetIso(0.5);

//   // Exemples : 3 points + 1 segment
//   implicit.AddPoint(Vector(-0.4,  0.0, 0.0), 0.35, 1.0, Color(0.9,0.2,0.2));
//   implicit.AddPoint(Vector( 0.4,  0.0, 0.0), 0.35, 1.0, Color(0.2,0.9,0.3));
//   implicit.AddPoint(Vector( 0.0,  0.4, 0.0), 0.35, 1.0, Color(0.2,0.6,0.9));
//   implicit.AddSegment(Vector(-0.6,-0.3, 0.0), Vector(0.6,-0.3,0.0),
//                       0.22, 1.0, Color(0.8,0.7,0.2));
// 	implicit.AddSegment(Vector(-0.6,-0.3, 0.0), Vector(0.0,0.6,0.0),
// 					  0.22, 1.0, Color(0.8,0.7,0.2));


//   Mesh implicitMesh;
//   implicit.Polygonize(64, implicitMesh, Box(2.0));

//   // Couleur simple (gris). Si tu veux colorer par primitive dominante,
//   // on peut garder le code de colorisation montré plus tôt.
//   std::vector<Color> cols(implicitMesh.Vertexes(), Color(0.8,0.8,0.8));

//   meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
//   UpdateGeometry();


    const double iso = 0.5;
    const double w = 1.0;
    const double R_values[] = {0.1, 0.2, 0.3, 0.4, 0.5};

    for (double R : R_values)
    {
        BlobField blob;
        blob.Clear();
        blob.SetIso(iso);

        // Trois points alignés
        blob.AddPoint(Vector(-0.4, 0.0, 0.0), R, w, Color(0.9, 0.2, 0.2));
        blob.AddPoint(Vector( 0.0, 0.0, 0.0), R, w, Color(0.2, 0.9, 0.3));
        blob.AddPoint(Vector( 0.4, 0.0, 0.0), R, w, Color(0.2, 0.3, 0.9));

        Mesh mesh;
        blob.Polygonize(128, mesh, Box(2.0));

        std::vector<Color> cols(mesh.Vertexes(), Color(0.8,0.8,0.8));
        meshColor = MeshColor(mesh, cols, mesh.VertexIndexes());
        UpdateGeometry();

        qDebug() << "R =" << R
                 << " => " << mesh.Vertexes() << " vertices,"
                 << mesh.Triangles() << " triangles.";
    }


    // const double iso = 0.5;
    // const double R = 0.3;

    // BlobField blob;
    // blob.Clear();
    // blob.SetIso(iso);

    // blob.AddPoint(Vector(-0.1, 0.0, 0.0), R, 1.0, Color(0.9, 0.2, 0.2)); // faible influence
    // blob.AddPoint(Vector( 0.1, 0.0, 0.0), R, 2.0, Color(0.2, 0.9, 0.3)); // poids plus fort

    // Mesh mesh;
    // blob.Polygonize(64, mesh, Box(2.0));

    // std::vector<Color> cols(mesh.Vertexes(), Color(0.8,0.8,0.8));
    // meshColor = MeshColor(mesh, cols, mesh.VertexIndexes());
    // UpdateGeometry();

    // qDebug() << "Weight test =>"
    //          << mesh.Vertexes() << "vertices,"
    //          << mesh.Triangles() << "triangles.";

}

void MainWindow::BenchmarkBlobs()
{
    const Box box(2.0);
    const int n = 128;
    std::vector<int> counts = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U(-0.8, 0.8);

    QString report = "Blobs Benchmark (n=" + QString::number(n) + ")\n";
    report += "count\tms\ttriangles\n";

    BlobField field;
    for (int m : counts)
    {
        
        field.SetIso(0.5); 

        for (int i = 0; i < m; ++i)
        {
            Vector p(U(rng), U(rng), U(rng));
            field.AddPoint(p, /*R*/0.25, /*w*/1.0,   // couleur = Color, pas Vector
                           Color(0.8, 0.8, 0.8));
        }

        Mesh mesh;
        QElapsedTimer t; t.start();
        field.Polygonize(n, mesh, box);
        qint64 ms = t.elapsed();

        report += QString("%1\t%2\t%3\n")
                    .arg(m)
                    .arg(ms)
                    .arg(mesh.Triangles());
    }

    // Après la boucle for :
    Mesh lastMesh;
    field.Polygonize(n, lastMesh, box);
    std::vector<Color> cols(lastMesh.Vertexes(), Color(0.8,0.8,0.8));
    meshColor = MeshColor(lastMesh, cols, lastMesh.VertexIndexes());
    UpdateGeometry();


    qDebug().noquote() << report;
    // QMessageBox::information(this, "Benchmark", report);
}

struct BenchRes { long long buildInc=0, buildBatch=0, evalInc=0, evalBatch=0; };

static long long time_ns(const std::function<void()>& fn, int repeat=1)
{
    long long best = LLONG_MAX;
    for (int i=0; i<repeat; ++i) {
        QElapsedTimer t; t.start();
        fn();
        long long ns = t.nsecsElapsed();
        if (ns < best) best = ns;
    }
    return best;
}

struct BenchImpact { Vector c; double R; double k; };

static BenchRes BenchErosion(int N=1000, int S=40, double k=0.02)
{
    std::vector<BenchImpact> impacts; impacts.reserve(N);
    for (int i=0; i<N; ++i) {
        double a = i * 0.123;
        impacts.push_back({ Vector(std::cos(a)*0.6, std::sin(a)*0.6, 0.0),
                            0.10 + ((i%7)==0 ? 0.04 : 0.0), k });
    }

    AnalyticScalarField base;

    ErosionIncremental inc(&base);
    long long nsBuildInc = time_ns([&]{
        for (const auto& s : impacts) inc.Apply(s.c, s.R, s.k);
    }, 1);

    ErosionField batch(&base);
    long long nsBuildBatch = time_ns([&]{
        for (const auto& s : impacts) batch.AddImpact(s.c, s.R, s.k);
    }, 1);

    auto Sample = [&](const AnalyticScalarField& f){
        volatile double acc=0.0;
        for (int z=0; z<S; ++z)
        for (int y=0; y<S; ++y)
        for (int x=0; x<S; ++x) {
            Vector p((x-S/2)/double(S/2), (y-S/2)/double(S/2), (z-S/2)/double(S/2));
            acc += f.Value(p);
        }
        return acc;
    };

    long long nsEvalInc = time_ns([&]{ (void)Sample(*inc.current); }, 3);
    long long nsEvalBatch = time_ns([&]{ (void)Sample(batch);          }, 3);

    BenchRes r;
    r.buildInc   = nsBuildInc   / 1000000;
    r.buildBatch = nsBuildBatch / 1000000;
    r.evalInc    = nsEvalInc    / 1000000;
    r.evalBatch  = nsEvalBatch  / 1000000;
    return r;
}

void MainWindow::RunBenchErosion()
{
    BenchRes R1 = BenchErosion(200,  40, 0.00);
    BenchRes R2 = BenchErosion(1000, 40, 0.02);
    BenchRes R3 = BenchErosion(3000, 40, 0.02);

    QString msg;
    auto line=[&](int N, const BenchRes& r){
        msg += QString("N=%1  Build: inc=%2 ms | batch=%3 ms   Eval: inc=%4 ms | batch=%5 ms\n")
               .arg(N).arg(r.buildInc).arg(r.buildBatch).arg(r.evalInc).arg(r.evalBatch);
    };
    line(200, R1); line(1000, R2); line(3000, R3);

    QMessageBox::information(this, "Bench Erosion", msg);
}

void MainWindow::ErosionImplicitExample()
{
    // PLUSIEURS CAS POUR LE RAPPORT DU TP QUE J'AI MIS DANS LE REPO ...

    // // BenchErosion();
    // // Base = une boule unité (AnalyticScalarField::Value = Norm(p)-1)
    // AnalyticScalarField base;

    // // Champ érodé
    // ErosionField eroded(&base);
    // eroded.AddImpact(Vector(0.3, 0.0, 0.0), 0.75, 0.08);
    // eroded.AddImpact(Vector(-0.2, 0.2, 0.0), 0.25, 0.00);

    // Mesh implicitMesh;
    // eroded.Polygonize(64, implicitMesh, Box(2.0)); // on réutilise TOUT ton pipeline Marching Cubes

    // std::vector<Color> cols(implicitMesh.Vertexes(), Color(0.8,0.8,0.8));
    // meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
    // UpdateGeometry();


    // // Base = sphère unité
    // AnalyticScalarField base;

    // // ErosionField = plusieurs sphères qui creusent
    // ErosionField eroded(&base);
    // eroded.AddImpact(Vector(0.85, 0.0, 0.0), 0.25, 0.08);
    // // eroded.AddImpact(Vector(0.0, -3.2, 0.0), 4.0, 0.08);
    // // eroded.AddImpact(Vector(-3.2, 0.0, 0.0), 4.0, 0.08);
    // // eroded.AddImpact(Vector(0.0, 3.2, 0.0), 4.0, 0.08); 
    // // eroded.AddImpact(Vector(0.0, 0.0, 3.2), 4.0, 0.08);
    // // eroded.AddImpact(Vector(0.0, 0.0, -3.2), 4.0, 0.08); 

    // Mesh m;
    // eroded.Polygonize(128, m, Box(2.0));

    // std::vector<Color> cols(m.Vertexes(), Color(0.8, 0.8, 0.8));
    // meshColor = MeshColor(m, cols, m.VertexIndexes());
    // UpdateGeometry();


    // AnalyticScalarField bas
    // ErosionField eroded(&base);

    // // anneau d’impacts qui se chevauchent => cavité qui perce
    // for (int i=0; i<10; ++i) {
    //     double a = (2.0*M_PI*i)/10.0;
    //     Vector c(0.85*std::cos(a), 0.85*std::sin(a), 0.0); // proche de la surface
    //     eroded.AddImpact(c, 0.3, 0.08, ImpactMode::Erode);
    // }

    // // un impact central plus petit pour "percer"
    // eroded.AddImpact(Vector(0.95, 0, 0), 0.22, 0.01, ImpactMode::Erode);

    // Mesh m; eroded.Polygonize(96, m, Box(2.5));
    // std::vector<Color> cols(m.Vertexes(), Color(0.8, 0.8, 0.8));
    // meshColor = MeshColor(m, cols, m.VertexIndexes());
    // UpdateGeometry();


    // AnalyticScalarField base;
    // ErosionField field(&base);
    // std::mt19937 rng(42);
    // std::uniform_real_distribution<double> U(-1.0, 1.0);
    // std::uniform_real_distribution<double> Rr(0.12, 0.35);

    // for (int i=0; i<40; ++i) {
    //     Vector dir = Normalized(Vector(U(rng), U(rng), U(rng)));
    //     double t   = 0.82 + 0.10*U(rng);        // autour de la surface
    //     Vector c   = dir * t;
    //     double R   = Rr(rng);
    //     double k   = 0.02 + 0.18*std::abs(U(rng));

    //     ImpactMode mode = (i%3==0) ? ImpactMode::Deposit : ImpactMode::Erode; // 2/3 creux, 1/3 bosses
    //     field.AddImpact(c, R, k, mode);
    // }

    // Mesh m; field.Polygonize(128, m, Box(2.5));
    // std::vector<Color> cols(m.Vertexes(), Color(0.8, 0.8, 0.8));
    // meshColor = MeshColor(m, cols, m.VertexIndexes());
    // UpdateGeometry();


    // AnalyticScalarField base;
    // ErosionField field(&base);

    // // coups nets (k=0) qui font des facettes/fractures
    // field.AddImpact(Vector(0.8, 0.0, 0.0), 0.45, 0.0, ImpactMode::Erode);
    // field.AddImpact(Vector(-0.7, 0.2, 0.0), 0.40, 0.0, ImpactMode::Erode);
    // field.AddImpact(Vector(0.0, -0.75, 0.2), 0.35, 0.0, ImpactMode::Erode);

    // // quelques dépôts durs
    // field.AddImpact(Vector(0.2, 0.75, -0.1), 0.28, 0.0, ImpactMode::Deposit);

    // Mesh m; field.Polygonize(96, m, Box(2.5));
    // std::vector<Color> cols(m.Vertexes(), Color(0.8, 0.8, 0.8));
    // meshColor = MeshColor(m, cols, m.VertexIndexes());
    // UpdateGeometry();


    struct BaseSphere : public AnalyticScalarField {
        double R0;
        explicit BaseSphere(double r) : R0(r) {}
        double Value(const Vector& p) const override { return Norm(p) - R0; }
    };

    const double R0 = 1.2;  
    BaseSphere base(R0);

    ErosionField eroded(&base);

    const double k_soft = 0.035;
    const double k_hard = 0.0;

    {
        const int N = 60;    
        const double rad = 0.14; 
        for (int i = 0; i < N; ++i)
        {
            double a = (2.0 * M_PI) * (i / double(N));
            Vector c(std::cos(a) * 0.85, 0.0, std::sin(a) * 0.85); 
            eroded.AddImpact(c, rad, k_soft); 
        }
    }

    {
        std::mt19937 rng(12345);
        std::uniform_real_distribution<double> Upos(-0.9, 0.9);
        std::uniform_real_distribution<double> Urad(0.08, 0.18);
        const int M = 200; // nombre d'impacts aléatoires
        for (int i = 0; i < M; ++i)
        {
            Vector c(Upos(rng), Upos(rng), Upos(rng));

            if (Norm(c) < 0.2) { i--; continue; }
            double R = Urad(rng);
            eroded.AddImpact(c, R, k_soft);
        }
    }

    {
        eroded.AddImpact(Vector( 0.95,  0.15,  0.10), 0.28, k_hard);
        eroded.AddImpact(Vector(-0.70, -0.25, -0.10), 0.32, k_hard);
        eroded.AddImpact(Vector( 0.10,  0.85, -0.20), 0.26, k_hard);
    }

    const int gridN = 128;     
    const double maxImpactR = 0.35;  
    const double extent = (R0 + maxImpactR) * 1.2;
    const Box box(extent);

    Mesh m;
    eroded.Polygonize(gridN, m, box);  

    std::vector<Color> cols(m.Vertexes(), Color(0.82, 0.82, 0.82));
    meshColor = MeshColor(m, cols, m.VertexIndexes());
    UpdateGeometry();

    qDebug() << "[Erosion v2] verts=" << m.Vertexes()
             << " tris=" << m.Triangles()
             << " gridN=" << gridN << " extent=" << extent;
}

static Mesh RotateMeshYX(const Mesh& m, double degY, double degX)
{
    const double ay = degY * M_PI / 180.0;
    const double ax = degX * M_PI / 180.0;
    const double cy = std::cos(ay), sy = std::sin(ay);
    const double cx = std::cos(ax), sx = std::sin(ax);

    std::vector<Vector> verts;
    verts.reserve(m.Vertexes());
    for (int i = 0; i < m.Vertexes(); ++i) {
        Vector p = m.Vertex(i);
        // rot Y
        double x =  cy*p[0] + sy*p[2];
        double z = -sy*p[0] + cy*p[2];
        // rot X
        double y =  cx*p[1] - sx*z;
        z        =  sx*p[1] + cx*z;
        verts.emplace_back(x, y, z);
    }

    std::vector<int> idx = m.VertexIndexes();

    Mesh r(verts, idx);
    r.SmoothNormals();
    return r;
}


void MainWindow::TerrainHeightmapExample()
{
    // PLUSIEURS CAS ...

    HeightParams P;

    // P.N = 193;        // taille raisonnable
    // P.scale = 0.25f;  // plan compact
    // P.amp = 40.0f;    // haut !
    // P.freq = 0.6f;    // (utilisé par fbm_ridged/detail)
    // P.octaves = 6;
    // P.lacunarity = 2.0f;
    // P.gain = 0.5f;
    // P.ridged = true;
    // P.warpStrength = 0.5f;
    // P.warpFreq = 0.4f;

    // P.N = 161;
    // P.scale = 0.22f;
    // P.amp = 55.0f;     // très haut
    // P.freq = 0.7f;
    // P.octaves = 6;
    // P.lacunarity = 2.1f;
    // P.gain = 0.48f;
    // P.ridged = true;
    // P.warpStrength = 0.7f;  // plis plus prononcés
    // P.warpFreq = 0.5f;

    // P.N = 193;
    // P.scale = 0.28f;
    // P.amp = 35.0f;
    // P.freq = 0.5f;
    // P.octaves = 5;
    // P.lacunarity = 1.9f;
    // P.gain = 0.55f;
    // P.ridged = false;   // collines plutôt que pics
    // P.warpStrength = 0.4f;
    // P.warpFreq = 0.35f;

    //     P.N = 193;
    // P.scale = 0.18f;
    // P.amp = 35.0f;
    // P.freq = 0.5f;
    // P.octaves = 5;
    // P.lacunarity = 1.9f;
    // P.gain = 0.55f;
    // P.ridged = false;   // collines plutôt que pics
    // P.warpStrength = 0.4f;
    // P.warpFreq = 0.35f;
    // P.heightBias = -5.0f; // Décalage vers le bas pour plus d'eau 


    // ——— TERRAIN SMOOTH, GROSSES COLLINES + GROS TROUS ———
    // P.N = 193;
    // P.scale = 0.25f;

    // P.amp   = 55.0f;    // ↑ plus de dynamique verticale
    // P.freq  = 0.45f;    // base pour fbm interne (détail modéré)
    // P.octaves = 5;
    // P.lacunarity = 2.0f;
    // P.gain  = 0.55f;

    // P.ridged = false;     // look smooth
    // P.warpStrength = 0.25f;
    // P.warpFreq     = 0.35f;

    // P.heightBias   = -6.0f; // ↓ décale tout vers le bas => zones < 0 (eau) garanties

    // P.N = 193;
    // P.scale = 0.22f;
    // P.amp = 35.0f;
    // P.freq = 0.45f;
    // P.octaves = 4;
    // P.lacunarity = 2.0f;
    // P.gain = 0.55f;
    // P.ridged = false;
    // P.warpStrength = 0.25f;
    // P.warpFreq = 0.35f;
    // P.heightBias = -4.0f;

    // P.model = NoiseType::Perlin;
    // P.useWarp   = true;  P.warpStrength = 0.3f;  P.warpFreq = 0.35f;
    // P.useBasins = true;  // grandes cuvettes (eau)
    // P.amp = 55.0f; P.freq = 0.45f; P.octaves = 5; P.lacunarity = 2.0f; P.gain = 0.52f;
    // P.heightBias = -6.0f;


    // P.N = 193;
    // P.scale = 0.25f;
    // P.amp = 50.0f;
    // P.freq = 0.5f;
    // P.octaves = 5;
    // P.lacunarity = 1.9f;
    // P.gain = 0.55f;
    // P.model = NoiseType::Perlin;
    // P.useWarp   = true;  P.warpStrength = 0.3f;  P.warpFreq = 0.3f;
    // P.useBasins = true;  // active les grandes zones “eau”
    // P.heightBias = -6.0f;


    P.N = 193;
    P.scale = 0.18f;
    P.amp = 35.0f;
    P.freq = 0.45f;
    P.octaves = 4;
    P.lacunarity = 2.1f;
    P.gain = 0.55f;
    P.model = NoiseType::Perlin;
    P.useWarp   = true;  P.warpStrength = 0.25f;  P.warpFreq = 0.35f;
    P.useBasins = false;
    P.heightBias = -2.0f;

    Mesh terrain;
    BuildTerrainMesh(P, terrain);

    // La j'ai annuler mais voila
    Mesh terrainRot = RotateMeshYX(terrain, /*yaw*/ 0.0, /*pitch*/ 0.0);

    const int n = terrainRot.Vertexes();
    std::vector<Color> cols(n);

    for (int i = 0; i < n; ++i)
    {
        const Vector& p = terrainRot.Vertex(i);
        float y = p[1];

        Color c;
        if (y < -10.0f)
            c = Color(0.05, 0.10, 0.40); // eau profonde
        else if (y < -9.0f)
            c = Color(0.10, 0.40, 0.60); // bord de mer
        else if (y < -8.0f)
            c = Color(0.10, 0.50, 0.10); // herbe
        else if (y < -7.0f)
            c = Color(0.50, 0.40, 0.30); // roche
        else
            c = Color(0.90, 0.90, 0.90); // neige

        cols[i] = c;
    }

    meshColor = MeshColor(terrainRot, cols, terrainRot.VertexIndexes());
    UpdateGeometry();
}



static void ReportMeshResidual(const SDFNode& node, const Mesh& m)
{
    const int n = m.Vertexes();
    double sum=0, sum2=0, mx=0;
    for (int i=0; i<n; ++i) {
        double r = std::abs(node.Value(m.Vertex(i)));
        sum += r; sum2 += r*r; if (r>mx) mx=r;
    }
    double mean = (n? sum/n : 0);
    double rms  = (n? std::sqrt(sum2/n) : 0);
    std::cout << "[Residual] mean="<<mean<<" rms="<<rms<<" max="<<mx<<"\n";
}

static void LogCaseHeader(const char* title, int gridN) {
    std::cout << "\n=== " << title << " (grid=" << gridN << ") ===\n";
}
static void RenderReport(MainWindow* self, const char* title, const SDFNode& node, int gridN=64)
{
    LogCaseHeader(title, gridN);
    Mesh m; node.Polygonize(gridN, m, Box(2.0));
    self->meshColor = MeshColor(m, std::vector<Color>(m.Vertexes(), Color(0.85,0.85,0.85)), m.VertexIndexes());
    self->UpdateGeometry();
    std::cout << "verts=" << m.Vertexes() << " tris=" << m.Triangles() << "\n";
    ReportMeshResidual(node, m);
}

#include <random>
#include <chrono>
#include <iostream>
using Clock = std::chrono::high_resolution_clock;

static void BenchSDF(const SDFNode& node, const char* label, int N = 10000000)
{
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(-1.5, 1.5);

    volatile double warm=0;
    for (int i=0; i<200000; ++i) {
        Vector p(uni(rng), uni(rng), uni(rng));
        warm += node.Value(p);
    }

    double acc = 0.0;
    auto t0 = Clock::now();
    for (int i=0; i<N; ++i) {
        Vector p(uni(rng), uni(rng), uni(rng));
        acc += node.Value(p);
    }
    auto t1 = Clock::now();

    double secs = std::chrono::duration<double>(t1-t0).count();
    double ns_per_call = (secs * 1e9) / double(N);
    double Mcalls_per_s = (double(N) / secs) / 1e6;

    std::cout.setf(std::ios::fixed); std::cout.precision(3);
    std::cout << "[Bench] " << label
              << " : " << secs << " s, "
              << ns_per_call << " ns/call, "
              << Mcalls_per_s << " M calls/s"
              << "  (checksum=" << acc << ")\n";
}

void MainWindow::SDFExample()
{
    using std::make_shared;

    auto sph = make_shared<SphereNode>(Vector(0,0,0), 0.7);
    auto box = make_shared<BoxNode>(Vector(0.35,0,0), Vector(0.35,0.28,0.25));
    auto cap = make_shared<CapsuleNode>(Vector(-0.5,0,0), Vector(0.5,0,0), 0.2);
    auto tor = make_shared<TorusNode>(Vector(0,0,0), 0.6, 0.18);
    auto uni = make_shared<UnionNode>(sph, box);             // min
    auto inter = make_shared<IntersectionNode>(sph, box);    // max
    auto diff = make_shared<DifferenceNode>(sph, box);       // max(a,-b)
    auto blend = make_shared<BlendNode>(sph, tor, 0.25);     // smooth union

    const int N = 10000000; // 10^7

    BenchSDF(*sph,    "Sphere",               N);
    BenchSDF(*box,    "Box",                  N);
    BenchSDF(*cap,    "Capsule",              N);
    BenchSDF(*tor,    "Torus",                N);
    BenchSDF(*uni,    "Union(Sphere,Box)",    N);
    BenchSDF(*inter,  "Intersection(S,B)",    N);
    BenchSDF(*diff,   "Difference(S-B)",      N);
    BenchSDF(*blend,  "Blend(Sphere,Torus)",  N);

    Mesh m; blend->Polygonize(64, m, Box(2.0));
    std::vector<Color> cols(m.Vertexes(), Color(0.8,0.8,0.8));
    meshColor = MeshColor(m, cols, m.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
