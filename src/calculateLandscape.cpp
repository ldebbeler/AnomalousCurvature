    std::vector<double> alphat(40);
    std::vector<double> apos(40);
    std::vector<double> aneg(40);

    for(std::size_t i=0; i<alphat.size(); i++){
        double at{ 4.0 - 0.05*i };
        newSe foo(at);
        alphat[i] = at;
        apos[i] = foo.Apos();
        aneg[i] = foo.Aneg();

        std::vector<double> x(201);
        std::vector<double> imag(201);
        std::vector<double> real(201);
        for(std::size_t i=0; i<x.size(); i++){
            double y{ -5.0 + i*0.05 };
            x[i] = y;
            imag[i] = foo.bubbleScale(y).imag();
            real[i] = foo.bubbleScale(y).real();
        }

        scalingValues bar;
        bar.m_x = x;
        bar.m_imag = imag;
        bar.m_real = real;

        writeNew write(bar);
        
        int ID = (int)( 400 - 5*i );
        std::string fileID{ std::to_string(ID) };

        std::string name{ prefix + fileID + suffix };

        H5::H5File file(name, H5F_ACC_TRUNC );
        writeNew writefile(bar);
        writefile.writeMainResults(file);
    }

    selfEnergyValues sev;
    sev.m_at = alphat;
    sev.m_apos = apos;
    sev.m_aneg = aneg;

    H5::H5File file(selfEnergy, H5F_ACC_TRUNC );
        writeNew writefile(sev);
        writefile.writeSe(file);

    std::cout << "File created: \t" << selfEnergy << '\n';
 
