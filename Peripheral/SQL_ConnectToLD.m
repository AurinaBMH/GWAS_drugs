function dbc = SQL_ConnectToLD()

% Open a database connection (hoping that the port forwarding has been set up
% locally to be able to use localhost:1234)
connSettings = struct();
connSettings.hostname = 'localhost';
connSettings.dbname = 'geneLD';
connSettings.username = 'root';
connSettings.password = 'ben1234';
dbc = SQL_opendatabase(connSettings);

end
