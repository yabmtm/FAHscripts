// target all code elements and add the example snippet to both
const snippets = document.querySelectorAll("pre code");

const code = `
04:23:07:I2:REQ127810:Client: 'rpvi8jm5bag9=0xd2297a3c5d285644@113.160.117.99' WORK_FAILED
04:23:07:W :REQ127810:Client reported Failed Assignment PRCG:14183(2,407,16)#22 from rpvi8jm5bag9=0xd2297a3c5d285644@113.160.117.99
04:23:07:I1:REQ127810:Retrying failed job PRCG:14183(2,407,16) (0 previous retries)
04:23:07:W :REQ127810:Empty WU results, nothing to save
04:23:07:I3:REQ127810:Response WORK_ACK
04:23:07:I1:Registering secure assignment key with CS at 128.252.203.2:443
04:23:07:I1:OUT127812:> POST https://128.252.203.2:443/ HTTP/1.1
04:23:07:I3:Connecting to 128.252.203.2:443
04:23:07:I1:OUT127812:< 128.252.203.2:443 HTTP/1.1 500 HTTP_INTERNAL_SERVER_ERROR
04:23:08:I1:REQ127806:< 98.26.1.60:53944 POST http://155.247.166.219/ HTTP/1.0
04:23:08:I2:REQ127806:Request WORK_RESULTS
04:23:08:I2:REQ127806:Client: 'siacovelli=0xa7ac28145d698d27@98.26.1.60' WORK_RESULTS
04:23:08:I1:REQ127806:Job PRCG:14196(19,111,64) accepted from siacovelli=0xa7ac28145d698d27@98.26.1.60
04:23:08:I3:REQ127806:Response WORK_ACK
04:23:08:I3:P14196:R19:C111:G65:Executing: /bin/sh -c "mv -f ../results64/checkpointState.xml ../state.xml"
04:23:09:I1:REQ127813:< 208.118.89.38:36781 POST http://155.247.166.219/ HTTP/1.0
04:23:09:I2:REQ127813:Request WORK_RESULTS
04:23:09:I2:REQ127813:Client: 've6jhc=0x457742a45db0e191@208.118.89.38' WORK_RESULTS
04:23:09:I1:REQ127813:Job PRCG:14188(1,536,11) accepted from ve6jhc=0x457742a45db0e191@208.118.89.38
04:23:09:I3:REQ127813:Response WORK_ACK
04:23:09:I3:P14188:R1:C536:G12:Executing: /bin/sh -c "/usr/local/bin/gromacs-5.0.4/bin/convert-tpr -s ../frame11.tpr -f ../results11/frame11.trr -o ../frame12.tpr -extend 2500"
04:23:10:I1:P14188:R1:C536:G12:GROMACS:    gmx convert-tpr, VERSION 5.0.4
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:GROMACS is written by:
04:23:10:I1:P14188:R1:C536:G12:Emile Apol         Rossen Apostolov   Herman J.C. Berendsen Par Bjelkmar       
04:23:10:I1:P14188:R1:C536:G12:Aldert van Buuren  Rudi van Drunen    Anton Feenstra     Sebastian Fritsch  
04:23:10:I1:P14188:R1:C536:G12:Gerrit Groenhof    Christoph Junghans Peter Kasson       Carsten Kutzner    
04:23:10:I1:P14188:R1:C536:G12:Per Larsson        Justin A. Lemkul   Magnus Lundborg    Pieter Meulenhoff  
04:23:10:I1:P14188:R1:C536:G12:Erik Marklund      Teemu Murtola      Szilard Pall       Sander Pronk       
04:23:10:I1:P14188:R1:C536:G12:Roland Schulz      Alexey Shvetsov    Michael Shirts     Alfons Sijbers     
04:23:10:I1:P14188:R1:C536:G12:Peter Tieleman     Christian Wennberg Maarten Wolf       
04:23:10:I1:P14188:R1:C536:G12:and the project leaders:
04:23:10:I1:P14188:R1:C536:G12:Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:Copyright (c) 1991-2000, University of Groningen, The Netherlands.
04:23:10:I1:P14188:R1:C536:G12:Copyright (c) 2001-2014, The GROMACS development team at
04:23:10:I1:P14188:R1:C536:G12:Uppsala University, Stockholm University and
04:23:10:I1:P14188:R1:C536:G12:the Royal Institute of Technology, Sweden.
04:23:10:I1:P14188:R1:C536:G12:check out http://www.gromacs.org for more information.
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:GROMACS is free software; you can redistribute it and/or modify it
04:23:10:I1:P14188:R1:C536:G12:under the terms of the GNU Lesser General Public License
04:23:10:I1:P14188:R1:C536:G12:as published by the Free Software Foundation; either version 2.1
04:23:10:I1:P14188:R1:C536:G12:of the License, or (at your option) any later version.
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:GROMACS:      gmx convert-tpr, VERSION 5.0.4
04:23:10:I1:P14188:R1:C536:G12:Executable:   /usr/local/bin/gromacs-5.0.4/bin/gmx
04:23:10:I1:P14188:R1:C536:G12:Library dir:  /usr/local/bin/gromacs-5.0.4/share/gromacs/top
04:23:10:I1:P14188:R1:C536:G12:Command line:
04:23:10:I1:P14188:R1:C536:G12:  convert-tpr -s ../frame11.tpr -f ../results11/frame11.trr -o ../frame12.tpr -extend 2500
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:Reading toplogy and stuff from ../frame11.tpr
04:23:10:I1:P14188:R1:C536:G12:Reading file ../frame11.tpr, VERSION 5.0.4 (single precision)
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:NOTE: Reading the state from trajectory is an obsolete feature of gmx convert-tpr.
04:23:10:I1:P14188:R1:C536:G12:      Continuation should be done by loading a checkpoint file with mdrun -cpi
04:23:10:I1:P14188:R1:C536:G12:      This guarantees that all state variables are transferred.
04:23:10:I1:P14188:R1:C536:G12:      gmx convert-tpr is now only useful for increasing nsteps,
04:23:10:I1:P14188:R1:C536:G12:      but even that can often be avoided by using mdrun -maxh
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:NOTE: The simulation uses pressure coupling and/or stochastic dynamics.
04:23:10:I1:P14188:R1:C536:G12:gmx convert-tpr can not provide binary identical continuation.
04:23:10:I1:P14188:R1:C536:G12:If you want that, supply a checkpoint file to mdrun
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:Changing ld-seed from 928811597 to 2198544166
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:READING COORDS, VELS AND BOX FROM TRAJECTORY ../results11/frame11.trr...
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:trn version: GMX_trn_file (single precision)
04:23:10:I1:P14188:R1:C536:G12:Read    trr frame      0: step 14000000 time  28000.000Read    trr frame      1: step 14500000 time  29000.000Read    trr frame      2: step 15000000 time  30000.000
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:Using frame of step 15000000 time 30000
04:23:10:I1:P14188:R1:C536:G12:Writing statusfile with starting step   15000000 and length    1250000 steps...
04:23:10:I1:P14188:R1:C536:G12:                                 time  30000.000 and length   2500.000 ps
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:gcq#377: "How will I know it's working right?" (MGMT)
04:23:10:I1:P14188:R1:C536:G12:
04:23:10:I1:P14188:R1:C536:G12:Extending remaining runtime of by 2500 ps (now 1250000 steps)
04:23:10:I1:REQ127815:< 185.97.199.106:36668 POST http://155.247.166.219/ HTTP/1.0
04:23:10:I2:REQ127815:Request WORK_FAILED
04:23:10:I2:REQ127815:Client: 'zniper=0x52dcdaac5d6837b2@185.97.199.106' WORK_FAILED
04:23:10:W :REQ127815:Client reported Failed Assignment PRCG:14196(5,336,74)#130 from zniper=0x52dcdaac5d6837b2@185.97.199.106
04:23:10:I1:REQ127815:Retrying failed job PRCG:14196(5,336,74) (1 previous retries)
04:23:10:W :REQ127815:Empty WU results, nothing to save
04:23:10:I3:REQ127815:Response WORK_ACK
04:23:11:I1:REQ127816:< 12.3.40.251:53294 POST http://155.247.166.219/ HTTP/1.0
04:23:11:I2:REQ127816:Request WORK_REQUEST
04:23:11:I2:REQ127816:Client: 'Anonymous=0x469348b85dd827b5@12.3.40.251' WORK_REQUEST
04:23:11:I1:REQ127816:Job PRCG:14190(13,87,33) assigned to Anonymous=0x469348b85dd827b5@12.3.40.251 with token 0x501aef5383588172
04:23:11:I3:REQ127816:Response WORK_ASSIGNMENT
04:23:11:I1:REQ127811:< 99.242.24.245:62414 POST http://155.247.166.219/ HTTP/1.0
04:23:11:I2:REQ127811:Request WORK_RESULTS
04:23:11:I2:REQ127811:Client: 'magakenny=0xe94f956c5889b96e@99.242.24.245' WORK_RESULTS
04:23:11:I1:REQ127811:Job PRCG:14182(8,174,15) accepted from magakenny=0xe94f956c5889b96e@99.242.24.245
04:23:11:I3:REQ127811:Response WORK_ACK
04:23:11:I3:P14182:R8:C174:G16:Executing: /bin/sh -c "/usr/local/bin/gromacs-5.0.4/bin/gmx convert-tpr -s ../frame15.tpr -f ../results15/frame15.trr      -o ../frame16.tpr -extend 5000"
04:23:12:I1:P14182:R8:C174:G16:GROMACS:    gmx convert-tpr, VERSION 5.0.4
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:GROMACS is written by:
04:23:12:I1:P14182:R8:C174:G16:Emile Apol         Rossen Apostolov   Herman J.C. Berendsen Par Bjelkmar       
04:23:12:I1:P14182:R8:C174:G16:Aldert van Buuren  Rudi van Drunen    Anton Feenstra     Sebastian Fritsch  
04:23:12:I1:P14182:R8:C174:G16:Gerrit Groenhof    Christoph Junghans Peter Kasson       Carsten Kutzner    
04:23:12:I1:P14182:R8:C174:G16:Per Larsson        Justin A. Lemkul   Magnus Lundborg    Pieter Meulenhoff  
04:23:12:I1:P14182:R8:C174:G16:Erik Marklund      Teemu Murtola      Szilard Pall       Sander Pronk       
04:23:12:I1:P14182:R8:C174:G16:Roland Schulz      Alexey Shvetsov    Michael Shirts     Alfons Sijbers     
04:23:12:I1:P14182:R8:C174:G16:Peter Tieleman     Christian Wennberg Maarten Wolf       
04:23:12:I1:P14182:R8:C174:G16:and the project leaders:
04:23:12:I1:P14182:R8:C174:G16:Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:Copyright (c) 1991-2000, University of Groningen, The Netherlands.
04:23:12:I1:P14182:R8:C174:G16:Copyright (c) 2001-2014, The GROMACS development team at
04:23:12:I1:P14182:R8:C174:G16:Uppsala University, Stockholm University and
04:23:12:I1:P14182:R8:C174:G16:the Royal Institute of Technology, Sweden.
04:23:12:I1:P14182:R8:C174:G16:check out http://www.gromacs.org for more information.
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:GROMACS is free software; you can redistribute it and/or modify it
04:23:12:I1:P14182:R8:C174:G16:under the terms of the GNU Lesser General Public License
04:23:12:I1:P14182:R8:C174:G16:as published by the Free Software Foundation; either version 2.1
04:23:12:I1:P14182:R8:C174:G16:of the License, or (at your option) any later version.
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:GROMACS:      gmx convert-tpr, VERSION 5.0.4
04:23:12:I1:P14182:R8:C174:G16:Executable:   /usr/local/bin/gromacs-5.0.4/bin/gmx
04:23:12:I1:P14182:R8:C174:G16:Library dir:  /usr/local/bin/gromacs-5.0.4/share/gromacs/top
04:23:12:I1:P14182:R8:C174:G16:Command line:
04:23:12:I1:P14182:R8:C174:G16:  gmx convert-tpr -s ../frame15.tpr -f ../results15/frame15.trr -o ../frame16.tpr -extend 5000
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:Reading toplogy and stuff from ../frame15.tpr
04:23:12:I1:P14182:R8:C174:G16:Reading file ../frame15.tpr, VERSION 5.0.4 (single precision)
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:NOTE: Reading the state from trajectory is an obsolete feature of gmx convert-tpr.
04:23:12:I1:P14182:R8:C174:G16:      Continuation should be done by loading a checkpoint file with mdrun -cpi
04:23:12:I1:P14182:R8:C174:G16:      This guarantees that all state variables are transferred.
04:23:12:I1:P14182:R8:C174:G16:      gmx convert-tpr is now only useful for increasing nsteps,
04:23:12:I1:P14182:R8:C174:G16:      but even that can often be avoided by using mdrun -maxh
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:READING COORDS, VELS AND BOX FROM TRAJECTORY ../results15/frame15.trr...
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:trn version: GMX_trn_file (single precision)
04:23:12:I1:P14182:R8:C174:G16:Read    trr frame      0: step 37500000 time  75000.000Read    trr frame      1: step 38000000 time  76000.000Read    trr frame      2: step 38500000 time  77000.000Read    trr frame      3: step 39000000 time  78000.000Read    trr frame      4: step 39500000 time  79000.000Read    trr frame      5: step 40000000 time  80000.000
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:Using frame of step 40000000 time 80000
04:23:12:I1:P14182:R8:C174:G16:Writing statusfile with starting step   40000000 and length    2500000 steps...
04:23:12:I1:P14182:R8:C174:G16:                                 time  80000.000 and length   5000.000 ps
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:gcq#62: "She Says She Can't Go Home Without a Chaperone" (E. Costello)
04:23:12:I1:P14182:R8:C174:G16:
04:23:12:I1:P14182:R8:C174:G16:Extending remaining runtime of by 5000 ps (now 2500000 steps)
04:23:14:I1:REQ127817:< 88.129.201.98:54376 POST http://155.247.166.219/ HTTP/1.0
04:23:14:I2:REQ127817:Request WORK_FAILED
04:23:14:I2:REQ127817:Client: 'biggels=0x9cc5252c5d757545@88.129.201.98' WORK_FAILED
04:23:14:W :REQ127817:Client reported Failed Assignment PRCG:14196(3,268,81)#156 from biggels=0x9cc5252c5d757545@88.129.201.98
04:23:14:I1:REQ127817:Retrying failed job PRCG:14196(3,268,81) (1 previous retries)
04:23:14:W :REQ127817:Empty WU results, nothing to save
04:23:14:I3:REQ127817:Response WORK_ACK
04:23:14:I1:REQ127818:< 199.212.27.247:56270 POST http://155.247.166.219/ HTTP/1.0
04:23:14:I2:REQ127818:Request WORK_FAILED
04:23:14:I2:REQ127818:Client: 'Ratman=0x7a11afec5da8cbe3@199.212.27.247' WORK_FAILED
04:23:14:W :REQ127818:Client reported Failed Assignment PRCG:14196(8,82,71)#134 from Ratman=0x7a11afec5da8cbe3@199.212.27.247
04:23:14:I1:REQ127818:Retrying failed job PRCG:14196(8,82,71) (2 previous retries)
04:23:14:W :REQ127818:Empty WU results, nothing to save
04:23:14:I3:REQ127818:Response WORK_ACK
04:23:24:I1:REQ127801:< 105.242.174.155:57387 POST http://155.247.166.219/ HTTP/1.0
04:23:24:I2:REQ127801:Request WORK_RESULTS
04:23:24:I2:REQ127801:Client: 'Bert_The_Dragon=0x5339f79228316762@105.242.174.155' WORK_RESULTS
04:23:25:I1:REQ127801:Job PRCG:14196(18,450,45) accepted from Bert_The_Dragon=0x5339f79228316762@105.242.174.155
04:23:25:I3:REQ127801:Response WORK_ACK
04:23:25:I3:P14196:R18:C450:G46:Executing: /bin/sh -c "mv -f ../results45/checkpointState.xml ../state.xml"
04:23:25:I1:REQ127820:< 109.252.80.202:15173 POST http://155.247.166.219/ HTTP/1.0
04:23:25:I2:REQ127820:Request WORK_REQUEST
04:23:25:I2:REQ127820:Client: 'Kostill_FLDC_1A2ad7gQ3ZcxDTpBRxrR7CP2L7Tzn8gtex=0xcc2008b85b41fe1b@109.252.80.202' WORK_REQUEST
04:23:25:I1:REQ127820:Job PRCG:14182(8,174,16) assigned to Kostill_FLDC_1A2ad7gQ3ZcxDTpBRxrR7CP2L7Tzn8gtex=0xcc2008b85b41fe1b@109.252.80.202 with token 0xead61638ad47f967
04:23:25:I3:REQ127820:Response WORK_ASSIGNMENT
04:23:27:I1:REQ127809:< 38.101.70.123:41822 POST http://155.247.166.219/ HTTP/1.0
04:23:27:I2:REQ127809:Request WORK_RESULTS
04:23:27:I2:REQ127809:Client: 'benjamindjb=0x482003605b9992a7@38.101.70.123' WORK_RESULTS
04:23:27:I1:REQ127809:Job PRCG:14183(1,150,15) accepted from benjamindjb=0x482003605b9992a7@38.101.70.123
04:23:27:I3:REQ127809:Response WORK_ACK
04:23:27:I3:P14183:R1:C150:G16:Executing: /bin/sh -c "/usr/local/bin/gromacs-5.0.4/bin/gmx convert-tpr -s ../frame15.tpr -f ../results15/frame15.trr      -o ../frame16.tpr -extend 5000"
04:23:27:I1:P14183:R1:C150:G16:GROMACS:    gmx convert-tpr, VERSION 5.0.4
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:GROMACS is written by:
04:23:27:I1:P14183:R1:C150:G16:Emile Apol         Rossen Apostolov   Herman J.C. Berendsen Par Bjelkmar       
04:23:27:I1:P14183:R1:C150:G16:Aldert van Buuren  Rudi van Drunen    Anton Feenstra     Sebastian Fritsch  
04:23:27:I1:P14183:R1:C150:G16:Gerrit Groenhof    Christoph Junghans Peter Kasson       Carsten Kutzner    
04:23:27:I1:P14183:R1:C150:G16:Per Larsson        Justin A. Lemkul   Magnus Lundborg    Pieter Meulenhoff  
04:23:27:I1:P14183:R1:C150:G16:Erik Marklund      Teemu Murtola      Szilard Pall       Sander Pronk       
04:23:27:I1:P14183:R1:C150:G16:Roland Schulz      Alexey Shvetsov    Michael Shirts     Alfons Sijbers     
04:23:27:I1:P14183:R1:C150:G16:Peter Tieleman     Christian Wennberg Maarten Wolf       
04:23:27:I1:P14183:R1:C150:G16:and the project leaders:
04:23:27:I1:P14183:R1:C150:G16:Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:Copyright (c) 1991-2000, University of Groningen, The Netherlands.
04:23:27:I1:P14183:R1:C150:G16:Copyright (c) 2001-2014, The GROMACS development team at
04:23:27:I1:P14183:R1:C150:G16:Uppsala University, Stockholm University and
04:23:27:I1:P14183:R1:C150:G16:the Royal Institute of Technology, Sweden.
04:23:27:I1:P14183:R1:C150:G16:check out http://www.gromacs.org for more information.
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:GROMACS is free software; you can redistribute it and/or modify it
04:23:27:I1:P14183:R1:C150:G16:under the terms of the GNU Lesser General Public License
04:23:27:I1:P14183:R1:C150:G16:as published by the Free Software Foundation; either version 2.1
04:23:27:I1:P14183:R1:C150:G16:of the License, or (at your option) any later version.
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:GROMACS:      gmx convert-tpr, VERSION 5.0.4
04:23:27:I1:P14183:R1:C150:G16:Executable:   /usr/local/bin/gromacs-5.0.4/bin/gmx
04:23:27:I1:P14183:R1:C150:G16:Library dir:  /usr/local/bin/gromacs-5.0.4/share/gromacs/top
04:23:27:I1:P14183:R1:C150:G16:Command line:
04:23:27:I1:P14183:R1:C150:G16:  gmx convert-tpr -s ../frame15.tpr -f ../results15/frame15.trr -o ../frame16.tpr -extend 5000
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:Reading toplogy and stuff from ../frame15.tpr
04:23:27:I1:P14183:R1:C150:G16:Reading file ../frame15.tpr, VERSION 5.0.4 (single precision)
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:NOTE: Reading the state from trajectory is an obsolete feature of gmx convert-tpr.
04:23:27:I1:P14183:R1:C150:G16:      Continuation should be done by loading a checkpoint file with mdrun -cpi
04:23:27:I1:P14183:R1:C150:G16:      This guarantees that all state variables are transferred.
04:23:27:I1:P14183:R1:C150:G16:      gmx convert-tpr is now only useful for increasing nsteps,
04:23:27:I1:P14183:R1:C150:G16:      but even that can often be avoided by using mdrun -maxh
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:READING COORDS, VELS AND BOX FROM TRAJECTORY ../results15/frame15.trr...
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:trn version: GMX_trn_file (single precision)
04:23:27:I1:P14183:R1:C150:G16:Read    trr frame      0: step 37500000 time  75000.000Read    trr frame      1: step 38000000 time  76000.000Read    trr frame      2: step 38500000 time  77000.000Read    trr frame      3: step 39000000 time  78000.000Read    trr frame      4: step 39500000 time  79000.000Read    trr frame      5: step 40000000 time  80000.000
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:Using frame of step 40000000 time 80000
04:23:27:I1:P14183:R1:C150:G16:Writing statusfile with starting step   40000000 and length    2500000 steps...
04:23:27:I1:P14183:R1:C150:G16:                                 time  80000.000 and length   5000.000 ps
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:gcq#295: "A Pretty Village Burning Makes a Pretty Fire" (David Sandstrom)
04:23:27:I1:P14183:R1:C150:G16:
04:23:27:I1:P14183:R1:C150:G16:Extending remaining runtime of by 5000 ps (now 2500000 steps)
04:23:34:I1:REQ127824:< 94.213.169.101:61674 POST http://155.247.166.219/ HTTP/1.0
04:23:34:I2:REQ127824:Request WORK_FAILED
04:23:34:I2:REQ127824:Client: 'turbodai=0x4a1b2f785dd14b9b@94.213.169.101' WORK_FAILED
04:23:34:W :REQ127824:Client reported Failed Assignment PRCG:14196(2,345,64)#111 from turbodai=0x4a1b2f785dd14b9b@94.213.169.101
04:23:34:I1:REQ127824:Retrying failed job PRCG:14196(2,345,64) (1 previous retries)
04:23:34:W :REQ127824:Empty WU results, nothing to save
04:23:34:I3:REQ127824:Response WORK_ACK
04:23:36:I1:REQ127825:< 72.198.216.172:50294 POST http://155.247.166.219/ HTTP/1.0
04:23:36:I2:REQ127825:Request WORK_FAILED
04:23:36:I2:REQ127825:Client: 'Maxim_Muir=0x583d95de5de47327@72.198.216.172' WORK_FAILED
04:23:36:W :REQ127825:Client reported Failed Assignment PRCG:14196(10,350,62)#113 from Maxim_Muir=0x583d95de5de47327@72.198.216.172
04:23:36:I1:REQ127825:Retrying failed job PRCG:14196(10,350,62) (0 previous retries)
04:23:36:W :REQ127825:Empty WU results, nothing to save
04:23:36:I3:REQ127825:Response WORK_ACK
04:23:39:I1:REQ127821:< 73.86.25.76:54854 POST http://155.247.166.219/ HTTP/1.0
04:23:39:I2:REQ127821:Request WORK_RESULTS
04:23:39:I2:REQ127821:Client: 'CoolGTX=0x350990e65d16ae83@73.86.25.76' WORK_RESULTS
04:23:40:I1:REQ127821:Job PRCG:14196(8,63,51) accepted from CoolGTX=0x350990e65d16ae83@73.86.25.76
04:23:40:I3:REQ127821:Response WORK_ACK
04:23:40:I3:P14196:R8:C63:G52:Executing: /bin/sh -c "mv -f ../results51/checkpointState.xml ../state.xml"
04:23:40:I1:REQ127715:< 176.58.225.149:40832 POST http://155.247.166.219/ HTTP/1.0
04:23:40:I2:REQ127715:Request WORK_RESULTS
04:23:40:I2:REQ127715:Client: 'Alex_Soilemezidis=0x352f60da5d3c1a42@176.58.225.149' WORK_RESULTS
04:23:41:I1:REQ127715:Job PRCG:14196(10,458,33) accepted from Alex_Soilemezidis=0x352f60da5d3c1a42@176.58.225.149
04:23:41:I3:REQ127715:Response WORK_ACK
04:23:41:I3:P14196:R10:C458:G34:Executing: /bin/sh -c "mv -f ../results33/checkpointState.xml ../state.xml"
04:23:42:I1:REQ127827:< 107.203.57.93:55486 POST http://155.247.166.219/ HTTP/1.0
04:23:42:I2:REQ127827:Request WORK_FAILED
04:23:42:I2:REQ127827:Client: 'patngayle=0xfe27765a5dac9960@107.203.57.93' WORK_FAILED
04:23:42:W :REQ127827:Client reported Failed Assignment PRCG:14196(0,173,63)#120 from patngayle=0xfe27765a5dac9960@107.203.57.93
04:23:42:I1:REQ127827:Retrying failed job PRCG:14196(0,173,63) (1 previous retries)
04:23:42:W :REQ127827:Empty WU results, nothing to save
04:23:42:I3:REQ127827:Response WORK_ACK
04:23:44:I1:REQ127814:< 199.116.115.135:50766 POST http://155.247.166.219/ HTTP/1.0
04:23:44:I2:REQ127814:Request WORK_RESULTS
04:23:44:I2:REQ127814:Client: 'AminItani=0x2e4945525b864859@199.116.115.135' WORK_RESULTS
04:23:44:I1:REQ127814:Job PRCG:14196(14,480,12) accepted from AminItani=0x2e4945525b864859@199.116.115.135
04:23:44:I3:REQ127814:Response WORK_ACK
04:23:44:I3:P14196:R14:C480:G13:Executing: /bin/sh -c "mv -f ../results12/checkpointState.xml ../state.xml"
04:23:44:I1:REQ127828:< 216.165.95.148:55815 POST http://155.247.166.219/ HTTP/1.0
04:23:44:I2:REQ127828:Request WORK_RESULTS
04:23:44:I2:REQ127828:Client: 'brizap=0xf417f2965a21c3db@216.165.95.148' WORK_RESULTS
04:23:45:I1:REQ127828:Job PRCG:14182(21,143,8) accepted from brizap=0xf417f2965a21c3db@216.165.95.148
04:23:45:I3:REQ127828:Response WORK_ACK
04:23:45:I3:P14182:R21:C143:G9:Executing: /bin/sh -c "/usr/local/bin/gromacs-5.0.4/bin/gmx convert-tpr -s ../frame8.tpr -f ../results8/frame8.trr      -o ../frame9.tpr -extend 5000"
04:23:45:I1:P14182:R21:C143:G9:GROMACS:    gmx convert-tpr, VERSION 5.0.4
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:GROMACS is written by:
04:23:45:I1:P14182:R21:C143:G9:Emile Apol         Rossen Apostolov   Herman J.C. Berendsen Par Bjelkmar       
04:23:45:I1:P14182:R21:C143:G9:Aldert van Buuren  Rudi van Drunen    Anton Feenstra     Sebastian Fritsch  
04:23:45:I1:P14182:R21:C143:G9:Gerrit Groenhof    Christoph Junghans Peter Kasson       Carsten Kutzner    
04:23:45:I1:P14182:R21:C143:G9:Per Larsson        Justin A. Lemkul   Magnus Lundborg    Pieter Meulenhoff  
04:23:45:I1:P14182:R21:C143:G9:Erik Marklund      Teemu Murtola      Szilard Pall       Sander Pronk       
04:23:45:I1:P14182:R21:C143:G9:Roland Schulz      Alexey Shvetsov    Michael Shirts     Alfons Sijbers     
04:23:45:I1:P14182:R21:C143:G9:Peter Tieleman     Christian Wennberg Maarten Wolf       
04:23:45:I1:P14182:R21:C143:G9:and the project leaders:
04:23:45:I1:P14182:R21:C143:G9:Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:Copyright (c) 1991-2000, University of Groningen, The Netherlands.
04:23:45:I1:P14182:R21:C143:G9:Copyright (c) 2001-2014, The GROMACS development team at
04:23:45:I1:P14182:R21:C143:G9:Uppsala University, Stockholm University and
04:23:45:I1:P14182:R21:C143:G9:the Royal Institute of Technology, Sweden.
04:23:45:I1:P14182:R21:C143:G9:check out http://www.gromacs.org for more information.
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:GROMACS is free software; you can redistribute it and/or modify it
04:23:45:I1:P14182:R21:C143:G9:under the terms of the GNU Lesser General Public License
04:23:45:I1:P14182:R21:C143:G9:as published by the Free Software Foundation; either version 2.1
04:23:45:I1:P14182:R21:C143:G9:of the License, or (at your option) any later version.
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:GROMACS:      gmx convert-tpr, VERSION 5.0.4
04:23:45:I1:P14182:R21:C143:G9:Executable:   /usr/local/bin/gromacs-5.0.4/bin/gmx
04:23:45:I1:P14182:R21:C143:G9:Library dir:  /usr/local/bin/gromacs-5.0.4/share/gromacs/top
04:23:45:I1:P14182:R21:C143:G9:Command line:
04:23:45:I1:P14182:R21:C143:G9:  gmx convert-tpr -s ../frame8.tpr -f ../results8/frame8.trr -o ../frame9.tpr -extend 5000
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:Reading toplogy and stuff from ../frame8.tpr
04:23:45:I1:P14182:R21:C143:G9:Reading file ../frame8.tpr, VERSION 5.0.4 (single precision)
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:NOTE: Reading the state from trajectory is an obsolete feature of gmx convert-tpr.
04:23:45:I1:P14182:R21:C143:G9:      Continuation should be done by loading a checkpoint file with mdrun -cpi
04:23:45:I1:P14182:R21:C143:G9:      This guarantees that all state variables are transferred.
04:23:45:I1:P14182:R21:C143:G9:      gmx convert-tpr is now only useful for increasing nsteps,
04:23:45:I1:P14182:R21:C143:G9:      but even that can often be avoided by using mdrun -maxh
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:READING COORDS, VELS AND BOX FROM TRAJECTORY ../results8/frame8.trr...
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:trn version: GMX_trn_file (single precision)
04:23:45:I1:P14182:R21:C143:G9:Read    trr frame      0: step 20000000 time  40000.000Read    trr frame      1: step 20500000 time  41000.000Read    trr frame      2: step 21000000 time  42000.000Read    trr frame      3: step 21500000 time  43000.000Read    trr frame      4: step 22000000 time  44000.000Read    trr frame      5: step 22500000 time  45000.000
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:Using frame of step 22500000 time 45000
04:23:45:I1:P14182:R21:C143:G9:Writing statusfile with starting step   22500000 and length    2500000 steps...
04:23:45:I1:P14182:R21:C143:G9:                                 time  45000.000 and length   5000.000 ps
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:gcq#190: "Ich Bin Ein Berliner" (J.F. Kennedy)
04:23:45:I1:P14182:R21:C143:G9:
04:23:45:I1:P14182:R21:C143:G9:Extending remaining runtime of by 5000 ps (now 2500000 steps)
04:23:45:I1:REQ127831:< 216.228.112.21:39173 POST http://155.247.166.219/ HTTP/1.0
04:23:45:I2:REQ127831:Request WORK_REQUEST
04:23:45:I2:REQ127831:Client: 'jtran00=0x0a2e0e465c5cce6f@216.228.112.21' WORK_REQUEST
04:23:45:I1:REQ127831:Job PRCG:14182(21,143,9) assigned to jtran00=0x0a2e0e465c5cce6f@216.228.112.21 with token 0xef0e90a40c657b5a
04:23:45:I3:REQ127831:Response WORK_ASSIGNMENT
04:23:45:I1:REQ127829:< 167.75.254.253:10400 POST http://155.247.166.219/ HTTP/1.0
04:23:45:I2:REQ127829:Request WORK_RESULTS
04:23:45:I2:REQ127829:Client: 'MungSu=0xc54d9d185a2e9b52@167.75.254.253' WORK_RESULTS
04:23:46:I1:REQ127829:Job PRCG:14190(10,135,14) accepted from MungSu=0xc54d9d185a2e9b52@167.75.254.253
04:23:46:I3:REQ127829:Response WORK_ACK
04:23:46:I3:P14190:R10:C135:G15:Executing: /bin/sh -c "/usr/local/bin/gromacs-5.0.4/bin/convert-tpr -s ../frame14.tpr -f ../results14/frame14.trr -o ../frame15.tpr -extend 2500"
04:23:46:I1:P14190:R10:C135:G15:GROMACS:    gmx convert-tpr, VERSION 5.0.4
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:GROMACS is written by:
04:23:46:I1:P14190:R10:C135:G15:Emile Apol         Rossen Apostolov   Herman J.C. Berendsen Par Bjelkmar       
04:23:46:I1:P14190:R10:C135:G15:Aldert van Buuren  Rudi van Drunen    Anton Feenstra     Sebastian Fritsch  
04:23:46:I1:P14190:R10:C135:G15:Gerrit Groenhof    Christoph Junghans Peter Kasson       Carsten Kutzner    
04:23:46:I1:P14190:R10:C135:G15:Per Larsson        Justin A. Lemkul   Magnus Lundborg    Pieter Meulenhoff  
04:23:46:I1:P14190:R10:C135:G15:Erik Marklund      Teemu Murtola      Szilard Pall       Sander Pronk       
04:23:46:I1:P14190:R10:C135:G15:Roland Schulz      Alexey Shvetsov    Michael Shirts     Alfons Sijbers     
04:23:46:I1:P14190:R10:C135:G15:Peter Tieleman     Christian Wennberg Maarten Wolf       
04:23:46:I1:P14190:R10:C135:G15:and the project leaders:
04:23:46:I1:P14190:R10:C135:G15:Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:Copyright (c) 1991-2000, University of Groningen, The Netherlands.
04:23:46:I1:P14190:R10:C135:G15:Copyright (c) 2001-2014, The GROMACS development team at
04:23:46:I1:P14190:R10:C135:G15:Uppsala University, Stockholm University and
04:23:46:I1:P14190:R10:C135:G15:the Royal Institute of Technology, Sweden.
04:23:46:I1:P14190:R10:C135:G15:check out http://www.gromacs.org for more information.
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:GROMACS is free software; you can redistribute it and/or modify it
04:23:46:I1:P14190:R10:C135:G15:under the terms of the GNU Lesser General Public License
04:23:46:I1:P14190:R10:C135:G15:as published by the Free Software Foundation; either version 2.1
04:23:46:I1:P14190:R10:C135:G15:of the License, or (at your option) any later version.
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:GROMACS:      gmx convert-tpr, VERSION 5.0.4
04:23:46:I1:P14190:R10:C135:G15:Executable:   /usr/local/bin/gromacs-5.0.4/bin/gmx
04:23:46:I1:P14190:R10:C135:G15:Library dir:  /usr/local/bin/gromacs-5.0.4/share/gromacs/top
04:23:46:I1:P14190:R10:C135:G15:Command line:
04:23:46:I1:P14190:R10:C135:G15:  convert-tpr -s ../frame14.tpr -f ../results14/frame14.trr -o ../frame15.tpr -extend 2500
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:Reading toplogy and stuff from ../frame14.tpr
04:23:46:I1:P14190:R10:C135:G15:Reading file ../frame14.tpr, VERSION 5.0.4 (single precision)
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:NOTE: Reading the state from trajectory is an obsolete feature of gmx convert-tpr.
04:23:46:I1:P14190:R10:C135:G15:      Continuation should be done by loading a checkpoint file with mdrun -cpi
04:23:46:I1:P14190:R10:C135:G15:      This guarantees that all state variables are transferred.
04:23:46:I1:P14190:R10:C135:G15:      gmx convert-tpr is now only useful for increasing nsteps,
04:23:46:I1:P14190:R10:C135:G15:      but even that can often be avoided by using mdrun -maxh
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:NOTE: The simulation uses pressure coupling and/or stochastic dynamics.
04:23:46:I1:P14190:R10:C135:G15:gmx convert-tpr can not provide binary identical continuation.
04:23:46:I1:P14190:R10:C135:G15:If you want that, supply a checkpoint file to mdrun
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:Changing ld-seed from 872307949 to 2654445794
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:READING COORDS, VELS AND BOX FROM TRAJECTORY ../results14/frame14.trr...
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:trn version: GMX_trn_file (single precision)
04:23:46:I1:P14190:R10:C135:G15:Read    trr frame      0: step 17500000 time  35000.000Read    trr frame      1: step 18000000 time  36000.000Read    trr frame      2: step 18500000 time  37000.000Read    trr frame      3: step 18750000 time  37500.000
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:Using frame of step 18750000 time 37500
04:23:46:I1:P14190:R10:C135:G15:Writing statusfile with starting step   18750000 and length    1250000 steps...
04:23:46:I1:P14190:R10:C135:G15:                                 time  37500.000 and length   2500.000 ps
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:gcq#156: "Way to Go Dude" (Beavis and Butthead)
04:23:46:I1:P14190:R10:C135:G15:
04:23:46:I1:P14190:R10:C135:G15:Extending remaining runtime of by 2500 ps (now 1250000 steps)
04:23:51:I1:REQ127834:< 76.104.26.20:55635 POST http://155.247.166.219/ HTTP/1.0
04:23:51:I2:REQ127834:Request WORK_FAILED
04:23:51:I2:REQ127834:Client: 'Ron_Michener=0xbfb62ac25bb8bbfc@76.104.26.20' WORK_FAILED
04:23:51:W :REQ127834:Client reported Failed Assignment PRCG:14196(16,267,58)#117 from Ron_Michener=0xbfb62ac25bb8bbfc@76.104.26.20
04:23:51:I1:REQ127834:Retrying failed job PRCG:14196(16,267,58) (2 previous retries)
04:23:51:W :REQ127834:Empty WU results, nothing to save
04:23:51:I3:REQ127834:Response WORK_ACK
04:23:52:I1:REQ127832:< 98.228.78.91:61180 POST http://155.247.166.219/ HTTP/1.0
04:23:52:I2:REQ127832:Request WORK_RESULTS
04:23:52:I2:REQ127832:Client: 'newtekie1=0x6a8f386c59b2efb0@98.228.78.91' WORK_RESULTS
04:23:52:I1:REQ127832:Job PRCG:14196(4,166,85) accepted from newtekie1=0x6a8f386c59b2efb0@98.228.78.91
04:23:52:I3:REQ127832:Response WORK_ACK
04:23:52:I3:P14196:R4:C166:G86:Executing: /bin/sh -c "mv -f ../results85/checkpointState.xml ../state.xml"
04:23:53:I1:REQ127835:< 82.203.172.213:40243 POST http://155.247.166.219/ HTTP/1.0
04:23:53:I2:REQ127835:Request WORK_FAILED
04:23:53:I2:REQ127835:Client: 'Kimmo_Laakia=0xcba52d0e5a2eb656@82.203.172.213' WORK_FAILED
04:23:53:W :REQ127835:Client reported Failed Assignment PRCG:14196(1,320,54)#92 from Kimmo_Laakia=0xcba52d0e5a2eb656@82.203.172.213
04:23:53:I1:REQ127835:Retrying failed job PRCG:14196(1,320,54) (0 previous retries)
04:23:53:W :REQ127835:Empty WU results, nothing to save
04:23:53:I3:REQ127835:Response WORK_ACK
04:23:56:I1:REQ127836:< 47.90.203.122:44982 POST http://155.247.166.219/ HTTP/1.0
04:23:56:I2:REQ127836:Request WORK_REQUEST
04:23:56:I2:REQ127836:Client: 'bnxf55aofrq4=0x99a4378e5d56e289@47.90.203.122' WORK_REQUEST
04:23:57:I1:REQ127836:Job PRCG:14182(5,92,4) assigned to bnxf55aofrq4=0x99a4378e5d56e289@47.90.203.122 with token 0x12e0b713b74745fb
04:23:57:I3:REQ127836:Response WORK_ASSIGNMENT
04:23:57:I1:REQ127837:< 199.116.115.135:64656 POST http://155.247.166.219/ HTTP/1.0
04:23:57:I2:REQ127837:Request WORK_FAILED
04:23:57:I2:REQ127837:Client: 'papanca=0x59e3da5e5dbb9134@199.116.115.135' WORK_FAILED
04:23:57:W :REQ127837:Client reported Failed Assignment PRCG:14196(12,389,94)#155 from papanca=0x59e3da5e5dbb9134@199.116.115.135
04:23:57:I1:REQ127837:Retrying failed job PRCG:14196(12,389,94) (0 previous retries)
04:23:57:W :REQ127837:Empty WU results, nothing to save
04:23:57:I3:REQ127837:Response WORK_ACK
04:23:58:I1:REQ127839:< 171.67.108.158:36247 POST https://155.247.166.219:443/ HTTP/1.1
04:23:58:I2:REQ127839:Secure Request REGISTER_KEY
04:23:58:W :REQ127839:REQ127839:171.67.108.158:36247:171.67.108.158 not allowed access to collection server
04:23:58:I1:REQ127838:< 195.171.160.231:39524 POST http://155.247.166.219/ HTTP/1.0
04:23:58:I2:REQ127838:Request WORK_RESULTS
04:23:58:I2:REQ127838:Client: 'MattL=0xf5aee6725ddffdf9@195.171.160.231' WORK_RESULTS
04:23:58:I1:REQ127838:Job PRCG:14182(7,99,13) accepted from MattL=0xf5aee6725ddffdf9@195.171.160.231
04:23:58:I3:REQ127838:Response WORK_ACK
04:23:58:I3:P14182:R7:C99:G14:Executing: /bin/sh -c "/usr/local/bin/gromacs-5.0.4/bin/gmx convert-tpr -s ../frame13.tpr -f ../results13/frame13.trr      -o ../frame14.tpr -extend 5000"
04:23:58:I1:P14182:R7:C99:G14:GROMACS:    gmx convert-tpr, VERSION 5.0.4
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:GROMACS is written by:
04:23:58:I1:P14182:R7:C99:G14:Emile Apol         Rossen Apostolov   Herman J.C. Berendsen Par Bjelkmar       
04:23:58:I1:P14182:R7:C99:G14:Aldert van Buuren  Rudi van Drunen    Anton Feenstra     Sebastian Fritsch  
04:23:58:I1:P14182:R7:C99:G14:Gerrit Groenhof    Christoph Junghans Peter Kasson       Carsten Kutzner    
04:23:58:I1:P14182:R7:C99:G14:Per Larsson        Justin A. Lemkul   Magnus Lundborg    Pieter Meulenhoff  
04:23:58:I1:P14182:R7:C99:G14:Erik Marklund      Teemu Murtola      Szilard Pall       Sander Pronk       
04:23:58:I1:P14182:R7:C99:G14:Roland Schulz      Alexey Shvetsov    Michael Shirts     Alfons Sijbers     
04:23:58:I1:P14182:R7:C99:G14:Peter Tieleman     Christian Wennberg Maarten Wolf       
04:23:58:I1:P14182:R7:C99:G14:and the project leaders:
04:23:58:I1:P14182:R7:C99:G14:Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:Copyright (c) 1991-2000, University of Groningen, The Netherlands.
04:23:58:I1:P14182:R7:C99:G14:Copyright (c) 2001-2014, The GROMACS development team at
04:23:58:I1:P14182:R7:C99:G14:Uppsala University, Stockholm University and
04:23:58:I1:P14182:R7:C99:G14:the Royal Institute of Technology, Sweden.
04:23:58:I1:P14182:R7:C99:G14:check out http://www.gromacs.org for more information.
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:GROMACS is free software; you can redistribute it and/or modify it
04:23:58:I1:P14182:R7:C99:G14:under the terms of the GNU Lesser General Public License
04:23:58:I1:P14182:R7:C99:G14:as published by the Free Software Foundation; either version 2.1
04:23:58:I1:P14182:R7:C99:G14:of the License, or (at your option) any later version.
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:GROMACS:      gmx convert-tpr, VERSION 5.0.4
04:23:58:I1:P14182:R7:C99:G14:Executable:   /usr/local/bin/gromacs-5.0.4/bin/gmx
04:23:58:I1:P14182:R7:C99:G14:Library dir:  /usr/local/bin/gromacs-5.0.4/share/gromacs/top
04:23:58:I1:P14182:R7:C99:G14:Command line:
04:23:58:I1:P14182:R7:C99:G14:  gmx convert-tpr -s ../frame13.tpr -f ../results13/frame13.trr -o ../frame14.tpr -extend 5000
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:Reading toplogy and stuff from ../frame13.tpr
04:23:58:I1:P14182:R7:C99:G14:Reading file ../frame13.tpr, VERSION 5.0.4 (single precision)
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:NOTE: Reading the state from trajectory is an obsolete feature of gmx convert-tpr.
04:23:58:I1:P14182:R7:C99:G14:      Continuation should be done by loading a checkpoint file with mdrun -cpi
04:23:58:I1:P14182:R7:C99:G14:      This guarantees that all state variables are transferred.
04:23:58:I1:P14182:R7:C99:G14:      gmx convert-tpr is now only useful for increasing nsteps,
04:23:58:I1:P14182:R7:C99:G14:      but even that can often be avoided by using mdrun -maxh
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:READING COORDS, VELS AND BOX FROM TRAJECTORY ../results13/frame13.trr...
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:trn version: GMX_trn_file (single precision)
04:23:58:I1:P14182:R7:C99:G14:Read    trr frame      0: step 32500000 time  65000.000Read    trr frame      1: step 33000000 time  66000.000Read    trr frame      2: step 33500000 time  67000.000Read    trr frame      3: step 34000000 time  68000.000Read    trr frame      4: step 34500000 time  69000.000Read    trr frame      5: step 35000000 time  70000.000
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:Using frame of step 35000000 time 70000
04:23:58:I1:P14182:R7:C99:G14:Writing statusfile with starting step   35000000 and length    2500000 steps...
04:23:58:I1:P14182:R7:C99:G14:                                 time  70000.000 and length   5000.000 ps
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:gcq#12: "Being Great is Not So Good" (Red Hot Chili Peppers)
04:23:58:I1:P14182:R7:C99:G14:
04:23:58:I1:P14182:R7:C99:G14:Extending remaining runtime of by 5000 ps (now 2500000 steps)
04:24:01:I1:REQ127842:< 2.37.163.197:55873 POST http://155.247.166.219/ HTTP/1.0
04:24:01:I2:REQ127842:Request WORK_FAILED
04:24:01:I2:REQ127842:Client: 'cafro=0xda1782ca5d09518a@2.37.163.197' WORK_FAILED
04:24:01:W :REQ127842:Client reported Failed Assignment PRCG:14196(1,403,51)#82 from cafro=0xda1782ca5d09518a@2.37.163.197
04:24:01:I1:REQ127842:Retrying failed job PRCG:14196(1,403,51) (0 previous retries)
04:24:01:W :REQ127842:Empty WU results, nothing to save
04:24:01:I3:REQ127842:Response WORK_ACK

`;

snippets.forEach(snippet => snippet.textContent = code);


// target the snippet in the regex section
const snippetRegex = document.querySelector(".snippet--regex pre code");

// include the regular expressions to find specific sections of text
// include them in an array of object, each detailing the expression and the class which needs to be applied upon the strings matching the connecting expression
const regularExpressions = [
        { 
          expression: /([0-9]+:[0-9]+:[0-9]+:..:|(REQ|OUT)[0-9]+:|^............#[0-9]+:)/gi,
          class: "timestamps"
        },
        { 
          expression: /([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+:[0-9]+.+|Connecting to|Client: '.+[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+'|(to |from ).+@.+|assign.+foldingathome.+|POST.+)/g,
          class: "ip_addrs"
        },
	{
          expression: /Executing.+/gi,
          class: "command"
	},
        { 
          expression: /(PRCG.[0-9]+.[0-9]+,[0-9]+,[0-9]+.|P[0-9]+.R[0-9]+.C[0-9]+.G[0-9]+|Job)/g,
          class: "projects"
        },
	{
	  expression: /Request.WORK_REQUEST|Request.WORK_RESULTS|WORK_REQUEST|WORK_RESULTS|WORK_ASSIGNMENT|WORK_ACK|Response|accepted|assigned|Updating.AS|KEY_ACCEPTED|Retrying f..... job/gi,
	  class: "good_response"
	},
	{
	  expression: /Request WORK_FAILED|WORK_FAILED|WORK_FAULTY|Client reported Failed Assignment|#[0-9]+|Request WORK_FAULTY|Core.+Assignment|gcq#.+/gi,
		class: "bad_response"
	},
	{
	  expression: /(Empty WU.+|Retrying failed job|.[0-9]+ previous retries\)|Secure.+|Registering.+CS.at|PLEASE_WAIT)/gi,
		class: "neutral_response"
	},
        {
          expression: /.+requesting shutdown|.+Clean exit/gi,
                class: "server_down"
        }
];

// loop through the array of regex and update the HTML structure of the snippet wrapping the strings matching the expressions in span elements, bearing a class paired to the expression itself
// class defined in the stylesheet to alter the appearance of the matching strings
regularExpressions.forEach((regularExpression) => snippetRegex.innerHTML = snippetRegex.innerHTML.replace(regularExpression.expression, `<span class=${regularExpression.class}>$&</span>`));
