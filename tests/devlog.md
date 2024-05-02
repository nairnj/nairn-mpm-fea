20244030 - Error en mandrel. La pieza no se estaba moviendo (material 3) y terminaba  en el filo del billet, 
           con lo cual no se le daba forma. Se le da cierta salida hacia afuera y se le impone velocidad como el pusher
         - El material 2 (para que el mandrel se desplace), no funciona ok si no se colocan también condiciones de borde 
           radiales.
         - Se agrega conductividad por contacto: https://osupdocs.forestry.oregonstate.edu/index.php/Thermal_Calculations#Conduction
         - LAS CONDICIONES DE BORDE TERMICAS NO FUNCIONAN CON TEMP EN EL BODY PARA CUERPOS RIGIDOS; Y SetTemperature CRASHEA
           ES MEJOR <TempBC adentro de BCLINE