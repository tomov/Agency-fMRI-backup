function goodSubjs = get_goodSubjs(which)

    switch (which)

        case 'S3'
            goodSubjs = [1 2 4 5 6 8 9 10 11 13 15 16 18 20 21 25 26 28 29 30 32 33 34  ];
            
        case 'S1'
            goodSubjs = [1 2 4 5 6 8 9 10 11 13 15 16 18 20 21 25 26 28 29 30 32 33 34       3    17    19    23    24    35    36    37 ];

        case 'S7'
            goodSubjs = [1 2 4 5 6 8 9 10 11 13 15 16 18 20 21 25 26 28 29 30 32 33 34       35 36];
            
        otherwise
            assert(false);

    end
